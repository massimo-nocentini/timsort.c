
/* 
	This is a glue c file for importing delta client c functions into Lua workflow.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <lua.h>
#include <lauxlib.h>
#include <time.h>



/* When we get into galloping mode, we stay there until both runs win less
 * often than MIN_GALLOP consecutive times.  See list_tsort.txt for more info.
 */
#define MIN_GALLOP 7

/* The maximum number of entries in a mergestate_t's pending-runs stack.
 * This is enough to sort arrays of size up to about
 *     32 * phi ** MAX_MERGE_PENDING
 * where phi ~= 1.618.  85 is ridiculouslylarge enough, good for an array
 * with 2**64 elements.
 */
#define MAX_MERGE_PENDING 85

/* Avoid malloc for small temp arrays. */
#define mergestate_t_TEMP_SIZE 256

typedef intptr_t signed_size_t;

typedef struct s_mergestate mergestate_t;

typedef int item_t;

typedef struct s_list_t {
    signed_size_t ob_size;
    
    /* Vector of pointers to list_t elements.  list_t[0] is ob_item[0], etc. */
    item_t *ob_item;

	lua_State *L;
    int table_absidx;
    
} list_t;


/* Lots of code for an adaptive, stable, natural mergesort.  There are many
 * pieces to this algorithm; read list_tsort.txt for overviews and details.
 */

/* A sortslice_t contains a pointer to an array of keys and a pointer to
 * an array of corresponding values.  In other words, keys[i]
 * corresponds with values[i].  If values == NULL, then the keys are
 * also the values.
 *
 * Several convenience routines are provided here, so that keys and
 * values are always moved in sync.
 */

typedef item_t * sortslice_t;

/* One mergestate_t exists on the stack per invocation of mergesort.  It's just
 * a convenient way to pass state around among the helper functions.
 */
struct s_slice {
    sortslice_t base;
    signed_size_t len;
};

struct s_mergestate {
    /* This controls when we get *into* galloping mode.  It's initialized
     * to MIN_GALLOP.  merge_lo and merge_hi tend to nudge it higher for
     * random data, and lower for highly structured data.
     */
    signed_size_t min_gallop;

    /* 'a' is temp storage to help with merges.  It contains room for
     * alloced entries.
     */
    sortslice_t a;        /* may point to temparray below */
    signed_size_t alloced;

    /* A stack of n pending runs yet to be merged.  Run #i starts at
     * address base[i] and extends for len[i] elements.  It's always
     * true (so long as the indices are in bounds) that
     *
     *     pending[i].base + pending[i].len == pending[i+1].base
     *
     * so we could cut the storage for this, but it's a minor amount,
     * and keeping all the info explicit simplifies the code.
     */
    int n;
    struct s_slice pending[MAX_MERGE_PENDING];

    /* 'a' points to this when possible, rather than muck with malloc. */
    item_t temparray[mergestate_t_TEMP_SIZE];

    list_t *list_tobject;
};

#define LOCAL_INLINE(type) static inline type

#define MERGE_GETMEM(MS, NEED) ((NEED) <= (MS)->alloced ? 0 : merge_getmem(MS, NEED))

/* Largest positive value of type signed_size_t. */
#define signed_size_t_MAX ((signed_size_t)(((size_t)-1)>>1))

/* Compare X to Y via "<".  Goto "fail" if the comparison raises an
   error.  Else "k" is set to true iff X<Y, and an "if (k)" block is
   started.  It makes more sense in context <wink>.  X and Y are item_t*s.
*/
#define IFLT(X, Y) if ((k = safe_object_compare(X, Y, ms)) < 0) goto fail; if (k)

LOCAL_INLINE(void) sortslice_t_memcpy(sortslice_t *s1, signed_size_t i, sortslice_t *s2, signed_size_t j, signed_size_t n)
{
    memcpy(*s1 + i, *s2 + j, sizeof(item_t ) * n);
}

LOCAL_INLINE(void) sortslice_t_copy_incr(sortslice_t *dst, sortslice_t *src)
{
    *(*dst)++ = *(*src)++;
}

LOCAL_INLINE(void) sortslice_t_advance(sortslice_t *slice, signed_size_t n)
{
    *slice += n;
}

LOCAL_INLINE(void) sortslice_t_memmove(sortslice_t *s1, signed_size_t i, sortslice_t *s2, signed_size_t j, signed_size_t n)
{
    memmove(*s1 + i, *s2 + j, sizeof(item_t ) * n);
}

LOCAL_INLINE(void) sortslice_t_copy(sortslice_t *s1, signed_size_t i, sortslice_t *s2, signed_size_t j)
{
    (*s1)[i] = (*s2)[j];
}

LOCAL_INLINE(void) sortslice_t_copy_decr(sortslice_t *dst, sortslice_t *src)
{
    *(*dst)-- = *(*src)--;
}

LOCAL_INLINE(int) safe_object_compare(item_t v, item_t w, mergestate_t *ms)
{
    int table_absidx = ms->list_tobject->table_absidx;

    lua_State *L = ms->list_tobject->L;
    lua_pushvalue(L, table_absidx + 2);
    lua_geti(L, table_absidx, v);
    lua_geti(L, table_absidx, w);

    lua_call(L, 2, 1);

    int lt = lua_toboolean(L, -1);

    lua_pop(L, 1);

    return lt;
}

/* Conceptually a mergestate_t's constructor. */
static void merge_init(mergestate_t *ms, signed_size_t list_t_size)
{
    assert(ms != NULL);
    ms->alloced = mergestate_t_TEMP_SIZE;
    ms->a  = ms->temparray;
    ms->n = 0;
    ms->min_gallop = MIN_GALLOP;
}

/* Reverse a slice of a list_t in place, from lo up to (exclusive) hi. */
static void reverse_slice(item_t *lo, item_t *hi)
{
    assert(lo && hi);

    --hi;
    while (lo < hi) {
        item_t t = *lo;
        *lo = *hi;
        *hi = t;
        ++lo;
        --hi;
    }
}

/* Compute a good value for the minimum run length; natural runs shorter
 * than this are boosted artificially via binary insertion.
 *
 * If n < 64, return n (it's too small to bother with fancy stuff).
 * Else if n is an exact power of 2, return 32.
 * Else return an int k, 32 <= k <= 64, such that n/k is close to, but
 * strictly less than, an exact power of 2.
 *
 * See list_tsort.txt for more info.
 */
static signed_size_t merge_compute_minrun(signed_size_t n)
{
    signed_size_t r = 0;           /* becomes 1 if any 1 bits are shifted off */

    assert(n >= 0);
    while (n >= 64) {
        r |= n & 1;
        n >>= 1;
    }
    return n + r;
}

/*
Return the length of the run beginning at lo, in the slice [lo, hi).  lo < hi
is required on entry.  "A run" is the longest ascending sequence, with

    lo[0] <= lo[1] <= lo[2] <= ...

or the longest descending sequence, with

    lo[0] > lo[1] > lo[2] > ...

Boolean *descending is set to 0 in the former case, or to 1 in the latter.
For its intended use in a stable mergesort, the strictness of the defn of
"descending" is needed so that the caller can safely reverse a descending
sequence without violating stability (strict > ensures there are no equal
elements to get out of order).

Returns -1 in case of error.
*/
static signed_size_t count_run(mergestate_t *ms, item_t *lo, item_t *hi, int *descending)
{
    signed_size_t k;
    signed_size_t n;

    assert(lo < hi);
    *descending = 0;
    ++lo;
    if (lo == hi)
        return 1;

    n = 2;
    IFLT(*lo, *(lo-1)) {
        *descending = 1;
        for (lo = lo+1; lo < hi; ++lo, ++n) {
            IFLT(*lo, *(lo-1))
                ;
            else
                break;
        }
    }
    else {
        for (lo = lo+1; lo < hi; ++lo, ++n) {
            IFLT(*lo, *(lo-1))
                break;
        }
    }

    return n;
fail:
    return -1;
}

static void reverse_sortslice_t(sortslice_t *s, signed_size_t n)
{
    reverse_slice(*s, *s + n);
}

/* binarysort is the best method for sorting small arrays: it does
   few compares, but can do data movement quadratic in the number of
   elements.
   [lo, hi) is a contiguous slice of a list_t, and is sorted via
   binary insertion.  This sort is stable.
   On entry, must have lo <= start <= hi, and that [lo, start) is already
   sorted (pass start == lo if you don't know!).
   If safe_object_compare(...) complains return -1, else 0.
   Even in case of error, the output slice will be some permutation of
   the input (nothing is lost or duplicated).
*/
static int binarysort(mergestate_t *ms, sortslice_t lo, item_t *hi, item_t *start)
{
    signed_size_t k;
    item_t *l, *p, *r;
    item_t pivot;

    assert(lo  <= start && start <= hi);
    /* assert [lo, start) is sorted */
    if (lo  == start)
        ++start;
    for (; start < hi; ++start) {
        /* set l to where *start belongs */
        l = lo ;
        r = start;
        pivot = *r;
        /* Invariants:
         * pivot >= all in [lo, l).
         * pivot  < all in [r, start).
         * The second is vacuously true at the start.
         */
        assert(l < r);
        do {
            p = l + ((r - l) >> 1);
            IFLT(pivot, *p)
                r = p;
            else
                l = p+1;
        } while (l < r);
        assert(l == r);
        /* The invariants still hold, so pivot >= all in [lo, l) and
           pivot < all in [l, start), so pivot belongs at l.  Note
           that if there are elements equal to pivot, l points to the
           first slot after them -- that's why this sort is stable.
           Slide over to make room.
           Caution: using memmove is much slower under MSVC 5;
           we're not usually moving many slots. */
        for (p = start; p > l; --p)
            *p = *(p-1);
        *l = pivot;
    }
    return 0;

 fail:
    return -1;
}


/*
Locate the proper position of key in a sorted vector; if the vector contains
an element equal to key, return the position immediately to the left of
the leftmost equal element.  [gallop_right() does the same except returns
the position to the right of the rightmost equal element (if any).]

"a" is a sorted vector with n elements, starting at a[0].  n must be > 0.

"hint" is an index at which to begin the search, 0 <= hint < n.  The closer
hint is to the final result, the faster this runs.

The return value is the int k in 0..n such that

    a[k-1] < key <= a[k]

pretending that *(a-1) is minus infinity and a[n] is plus infinity.  IOW,
key belongs at index k; or, IOW, the first k elements of a should precede
key, and the last n-k should follow key.

Returns -1 on error.  See list_tsort.txt for info on the method.
*/
static signed_size_t gallop_left(mergestate_t *ms, item_t key, item_t *a, signed_size_t n, signed_size_t hint)
{
    signed_size_t ofs;
    signed_size_t lastofs;
    signed_size_t k;

    assert(key && a && n > 0 && hint >= 0 && hint < n);

    a += hint;
    lastofs = 0;
    ofs = 1;
    IFLT(*a, key) {
        /* a[hint] < key -- gallop right, until
         * a[hint + lastofs] < key <= a[hint + ofs]
         */
        const signed_size_t maxofs = n - hint;             /* &a[n-1] is highest */
        while (ofs < maxofs) {
            IFLT(a[ofs], key) {
                lastofs = ofs;
                assert(ofs <= (signed_size_t_MAX - 1) / 2);
                ofs = (ofs << 1) + 1;
            }
            else                /* key <= a[hint + ofs] */
                break;
        }
        if (ofs > maxofs)
            ofs = maxofs;
        /* Translate back to offsets relative to &a[0]. */
        lastofs += hint;
        ofs += hint;
    }
    else {
        /* key <= a[hint] -- gallop left, until
         * a[hint - ofs] < key <= a[hint - lastofs]
         */
        const signed_size_t maxofs = hint + 1;             /* &a[0] is lowest */
        while (ofs < maxofs) {
            IFLT(*(a-ofs), key)
                break;
            /* key <= a[hint - ofs] */
            lastofs = ofs;
            assert(ofs <= (signed_size_t_MAX - 1) / 2);
            ofs = (ofs << 1) + 1;
        }
        if (ofs > maxofs)
            ofs = maxofs;
        /* Translate back to positive offsets relative to &a[0]. */
        k = lastofs;
        lastofs = hint - ofs;
        ofs = hint - k;
    }
    a -= hint;

    assert(-1 <= lastofs && lastofs < ofs && ofs <= n);
    /* Now a[lastofs] < key <= a[ofs], so key belongs somewhere to the
     * right of lastofs but no farther right than ofs.  Do a binary
     * search, with invariant a[lastofs-1] < key <= a[ofs].
     */
    ++lastofs;
    while (lastofs < ofs) {
        signed_size_t m = lastofs + ((ofs - lastofs) >> 1);

        IFLT(a[m], key)
            lastofs = m+1;              /* a[m] < key */
        else
            ofs = m;                    /* key <= a[m] */
    }
    assert(lastofs == ofs);             /* so a[ofs-1] < key <= a[ofs] */
    return ofs;

fail:
    return -1;
}

/*
Exactly like gallop_left(), except that if key already exists in a[0:n],
finds the position immediately to the right of the rightmost equal value.

The return value is the int k in 0..n such that

    a[k-1] <= key < a[k]

or -1 if error.

The code duplication is massive, but this is enough different given that
we're sticking to "<" comparisons that it's much harder to follow if
written as one routine with yet another "left or right?" flag.
*/
static signed_size_t gallop_right (mergestate_t *ms, item_t key, item_t *a, signed_size_t n, signed_size_t hint)
{
    signed_size_t ofs;
    signed_size_t lastofs;
    signed_size_t k;

    assert(key && a && n > 0 && hint >= 0 && hint < n);

    a += hint;
    lastofs = 0;
    ofs = 1;
    IFLT(key, *a) {
        /* key < a[hint] -- gallop left, until
         * a[hint - ofs] <= key < a[hint - lastofs]
         */
        const signed_size_t maxofs = hint + 1;             /* &a[0] is lowest */
        while (ofs < maxofs) {
            IFLT(key, *(a-ofs)) {
                lastofs = ofs;
                assert(ofs <= (signed_size_t_MAX - 1) / 2);
                ofs = (ofs << 1) + 1;
            }
            else                /* a[hint - ofs] <= key */
                break;
        }
        if (ofs > maxofs)
            ofs = maxofs;
        /* Translate back to positive offsets relative to &a[0]. */
        k = lastofs;
        lastofs = hint - ofs;
        ofs = hint - k;
    }
    else {
        /* a[hint] <= key -- gallop right, until
         * a[hint + lastofs] <= key < a[hint + ofs]
        */
        const signed_size_t maxofs = n - hint;             /* &a[n-1] is highest */
        while (ofs < maxofs) {
            IFLT(key, a[ofs])
                break;
            /* a[hint + ofs] <= key */
            lastofs = ofs;
            assert(ofs <= (signed_size_t_MAX - 1) / 2);
            ofs = (ofs << 1) + 1;
        }
        if (ofs > maxofs)
            ofs = maxofs;
        /* Translate back to offsets relative to &a[0]. */
        lastofs += hint;
        ofs += hint;
    }
    a -= hint;

    assert(-1 <= lastofs && lastofs < ofs && ofs <= n);
    /* Now a[lastofs] <= key < a[ofs], so key belongs somewhere to the
     * right of lastofs but no farther right than ofs.  Do a binary
     * search, with invariant a[lastofs-1] <= key < a[ofs].
     */
    ++lastofs;
    while (lastofs < ofs) {
        signed_size_t m = lastofs + ((ofs - lastofs) >> 1);

        IFLT(key, a[m])
            ofs = m;                    /* key < a[m] */
        else
            lastofs = m+1;              /* a[m] <= key */
    }
    assert(lastofs == ofs);             /* so a[ofs-1] <= key < a[ofs] */
    return ofs;

fail:
    return -1;
}

/* Free all the temp memory owned by the mergestate_t.  This must be called
 * when you're done with a mergestate_t, and may be called before then if
 * you want to free the temp memory early.
 */
static void merge_freemem(mergestate_t *ms)
{
    assert(ms != NULL);
    if (ms->a  != ms->temparray) {
        free(ms->a );
        ms->a  = NULL;
    }
}

void * checked_malloc(size_t size)
{
    /* see PyMem_RawMalloc() */
    if (size > (size_t)signed_size_t_MAX)
        return NULL;
    //return _PyMem.malloc(_PyMem.ctx, size);
    return malloc(size);
}

/* Ensure enough temp memory for 'need' array slots is available.
 * Returns 0 on success and -1 if the memory can't be gotten.
 */
static int merge_getmem(mergestate_t *ms, signed_size_t need)
{
    assert(ms != NULL);
    if (need <= ms->alloced)
        return 0;
    
    /* Don't realloc!  That can cost cycles to copy the old data, but
     * we don't care what's in the block.
     */
    merge_freemem(ms);
    if ((size_t)need > signed_size_t_MAX / sizeof(item_t )) {
        //PyErr_NoMemory();
        return -1;
    }
    ms->a  = (item_t *)checked_malloc(need * sizeof(item_t ));
    if (ms->a  != NULL) {
        ms->alloced = need;
        
        return 0;
    }
    //PyErr_NoMemory();
    return -1;
}

/* Merge the na elements starting at ssa with the nb elements starting at
 * ssb  = ssa  + na in a stable way, in-place.  na and nb must be > 0.
 * Must also have that ssa [na-1] belongs at the end of the merge, and
 * should have na <= nb.  See list_tsort.txt for more info.  Return 0 if
 * successful, -1 if error.
 */
static signed_size_t merge_lo(mergestate_t *ms, sortslice_t ssa, signed_size_t na, sortslice_t ssb, signed_size_t nb)
{
    signed_size_t k;
    sortslice_t dest;
    int result = -1;            /* guilty until proved innocent */
    signed_size_t min_gallop;

    assert(ms && ssa  && ssb  && na > 0 && nb > 0);
    assert(ssa  + na == ssb );
    if (MERGE_GETMEM(ms, na) < 0)
        return -1;
    sortslice_t_memcpy(&ms->a, 0, &ssa, 0, na);
    dest = ssa;
    ssa = ms->a;

    sortslice_t_copy_incr(&dest, &ssb);
    --nb;
    if (nb == 0)
        goto Succeed;
    if (na == 1)
        goto CopyB;

    min_gallop = ms->min_gallop;
    for (;;) {
        signed_size_t acount = 0;          /* # of times A won in a row */
        signed_size_t bcount = 0;          /* # of times B won in a row */

        /* Do the straightforward thing until (if ever) one run
         * appears to win consistently.
         */
        for (;;) {
            assert(na > 1 && nb > 0);
            k = safe_object_compare(ssb [0], ssa [0], ms);
            if (k) {
                if (k < 0)
                    goto Fail;
                sortslice_t_copy_incr(&dest, &ssb);
                ++bcount;
                acount = 0;
                --nb;
                if (nb == 0)
                    goto Succeed;
                if (bcount >= min_gallop)
                    break;
            }
            else {
                sortslice_t_copy_incr(&dest, &ssa);
                ++acount;
                bcount = 0;
                --na;
                if (na == 1)
                    goto CopyB;
                if (acount >= min_gallop)
                    break;
            }
        }

        /* One run is winning so consistently that galloping may
         * be a huge win.  So try that, and continue galloping until
         * (if ever) neither run appears to be winning consistently
         * anymore.
         */
        ++min_gallop;
        do {
            assert(na > 1 && nb > 0);
            min_gallop -= min_gallop > 1;
            ms->min_gallop = min_gallop;
            k = gallop_right(ms, ssb [0], ssa , na, 0);
            acount = k;
            if (k) {
                if (k < 0)
                    goto Fail;
                sortslice_t_memcpy(&dest, 0, &ssa, 0, k);
                sortslice_t_advance(&dest, k);
                sortslice_t_advance(&ssa, k);
                na -= k;
                if (na == 1)
                    goto CopyB;
                /* na==0 is impossible now if the comparison
                 * function is consistent, but we can't assume
                 * that it is.
                 */
                if (na == 0)
                    goto Succeed;
            }
            sortslice_t_copy_incr(&dest, &ssb);
            --nb;
            if (nb == 0)
                goto Succeed;

            k = gallop_left(ms, ssa [0], ssb , nb, 0);
            bcount = k;
            if (k) {
                if (k < 0)
                    goto Fail;
                sortslice_t_memmove(&dest, 0, &ssb, 0, k);
                sortslice_t_advance(&dest, k);
                sortslice_t_advance(&ssb, k);
                nb -= k;
                if (nb == 0)
                    goto Succeed;
            }
            sortslice_t_copy_incr(&dest, &ssa);
            --na;
            if (na == 1)
                goto CopyB;
        } while (acount >= MIN_GALLOP || bcount >= MIN_GALLOP);
        ++min_gallop;           /* penalize it for leaving galloping mode */
        ms->min_gallop = min_gallop;
    }
Succeed:
    result = 0;
Fail:
    if (na)
        sortslice_t_memcpy(&dest, 0, &ssa, 0, na);
    return result;
CopyB:
    assert(na == 1 && nb > 0);
    /* The last element of ssa belongs at the end of the merge. */
    sortslice_t_memmove(&dest, 0, &ssb, 0, nb);
    sortslice_t_copy(&dest, nb, &ssa, 0);
    return 0;
}


/* Merge the na elements starting at pa with the nb elements starting at
 * ssb  = ssa  + na in a stable way, in-place.  na and nb must be > 0.
 * Must also have that ssa [na-1] belongs at the end of the merge, and
 * should have na >= nb.  See list_tsort.txt for more info.  Return 0 if
 * successful, -1 if error.
 */
static signed_size_t merge_hi(mergestate_t *ms, sortslice_t ssa, signed_size_t na, sortslice_t ssb, signed_size_t nb)
{
    signed_size_t k;
    sortslice_t dest, basea, baseb;
    int result = -1;            /* guilty until proved innocent */
    signed_size_t min_gallop;

    assert(ms && ssa  && ssb  && na > 0 && nb > 0);
    assert(ssa  + na == ssb );
    if (MERGE_GETMEM(ms, nb) < 0)
        return -1;
    dest = ssb;
    sortslice_t_advance(&dest, nb-1);
    sortslice_t_memcpy(&ms->a, 0, &ssb, 0, nb);
    basea = ssa;
    baseb = ms->a;
    ssb  = ms->a  + nb - 1;
    
    sortslice_t_advance(&ssa, na - 1);

    sortslice_t_copy_decr(&dest, &ssa);
    --na;
    if (na == 0)
        goto Succeed;
    if (nb == 1)
        goto CopyA;

    min_gallop = ms->min_gallop;
    for (;;) {
        signed_size_t acount = 0;          /* # of times A won in a row */
        signed_size_t bcount = 0;          /* # of times B won in a row */

        /* Do the straightforward thing until (if ever) one run
         * appears to win consistently.
         */
        for (;;) {
            assert(na > 0 && nb > 1);
            k = safe_object_compare(ssb [0], ssa [0], ms);
            if (k) {
                if (k < 0)
                    goto Fail;
                sortslice_t_copy_decr(&dest, &ssa);
                ++acount;
                bcount = 0;
                --na;
                if (na == 0)
                    goto Succeed;
                if (acount >= min_gallop)
                    break;
            }
            else {
                sortslice_t_copy_decr(&dest, &ssb);
                ++bcount;
                acount = 0;
                --nb;
                if (nb == 1)
                    goto CopyA;
                if (bcount >= min_gallop)
                    break;
            }
        }

        /* One run is winning so consistently that galloping may
         * be a huge win.  So try that, and continue galloping until
         * (if ever) neither run appears to be winning consistently
         * anymore.
         */
        ++min_gallop;
        do {
            assert(na > 0 && nb > 1);
            min_gallop -= min_gallop > 1;
            ms->min_gallop = min_gallop;
            k = gallop_right(ms, ssb [0], basea , na, na-1);
            if (k < 0)
                goto Fail;
            k = na - k;
            acount = k;
            if (k) {
                sortslice_t_advance(&dest, -k);
                sortslice_t_advance(&ssa, -k);
                sortslice_t_memmove(&dest, 1, &ssa, 1, k);
                na -= k;
                if (na == 0)
                    goto Succeed;
            }
            sortslice_t_copy_decr(&dest, &ssb);
            --nb;
            if (nb == 1)
                goto CopyA;

            k = gallop_left(ms, ssa [0], baseb , nb, nb-1);
            if (k < 0)
                goto Fail;
            k = nb - k;
            bcount = k;
            if (k) {
                sortslice_t_advance(&dest, -k);
                sortslice_t_advance(&ssb, -k);
                sortslice_t_memcpy(&dest, 1, &ssb, 1, k);
                nb -= k;
                if (nb == 1)
                    goto CopyA;
                /* nb==0 is impossible now if the comparison
                 * function is consistent, but we can't assume
                 * that it is.
                 */
                if (nb == 0)
                    goto Succeed;
            }
            sortslice_t_copy_decr(&dest, &ssa);
            --na;
            if (na == 0)
                goto Succeed;
        } while (acount >= MIN_GALLOP || bcount >= MIN_GALLOP);
        ++min_gallop;           /* penalize it for leaving galloping mode */
        ms->min_gallop = min_gallop;
    }
Succeed:
    result = 0;
Fail:
    if (nb)
        sortslice_t_memcpy(&dest, -(nb-1), &baseb, 0, nb);
    return result;
CopyA:
    assert(nb == 1 && na > 0);
    /* The first element of ssb belongs at the front of the merge. */
    sortslice_t_memmove(&dest, 1-na, &ssa, 1-na, na);
    sortslice_t_advance(&dest, -na);
    sortslice_t_advance(&ssa, -na);
    sortslice_t_copy(&dest, 0, &ssb, 0);
    return 0;
}

/* Merge the two runs at stack indices i and i+1.
 * Returns 0 on success, -1 on error.
 */
static signed_size_t merge_at (mergestate_t *ms, signed_size_t i)
{
    sortslice_t ssa, ssb;
    signed_size_t na, nb;
    signed_size_t k;

    assert(ms != NULL);
    assert(ms->n >= 2);
    assert(i >= 0);
    assert(i == ms->n - 2 || i == ms->n - 3);

    ssa = ms->pending[i].base;
    na = ms->pending[i].len;
    ssb = ms->pending[i+1].base;
    nb = ms->pending[i+1].len;
    assert(na > 0 && nb > 0);
    assert(ssa  + na == ssb );

    /* Record the length of the combined runs; if i is the 3rd-last
     * run now, also slide over the last run (which isn't involved
     * in this merge).  The current run i+1 goes away in any case.
     */
    ms->pending[i].len = na + nb;
    if (i == ms->n - 3)
        ms->pending[i+1] = ms->pending[i+2];
    --ms->n;

    /* Where does b start in a?  Elements in a before that can be
     * ignored (already in place).
     */
    k = gallop_right(ms, *ssb , ssa , na, 0);
    if (k < 0)
        return -1;
    sortslice_t_advance(&ssa, k);
    na -= k;
    if (na == 0)
        return 0;

    /* Where does a end in b?  Elements in b after that can be
     * ignored (already in place).
     */
    nb = gallop_left(ms, ssa [na-1], ssb , nb, nb-1);
    if (nb <= 0)
        return nb;

    /* Merge what remains of the runs, using a temp array with
     * min(na, nb) elements.
     */
    if (na <= nb)
        return merge_lo(ms, ssa, na, ssb, nb);
    else
        return merge_hi(ms, ssa, na, ssb, nb);
}


/* Regardless of invariants, merge all runs on the stack until only one
 * remains.  This is used at the end of the mergesort.
 *
 * Returns 0 on success, -1 on error.
 */
static int merge_force_collapse(mergestate_t *ms)
{
    struct s_slice *p = ms->pending;

    assert(ms);
    while (ms->n > 1) {
        signed_size_t n = ms->n - 2;
        if (n > 0 && p[n-1].len < p[n+1].len)
            --n;
        if (merge_at(ms, n) < 0)
            return -1;
    }
    return 0;
}

/* Examine the stack of runs waiting to be merged, merging adjacent runs
 * until the stack invariants are re-established:
 *
 * 1. len[-3] > len[-2] + len[-1]
 * 2. len[-2] > len[-1]
 *
 * See list_tsort.txt for more info.
 *
 * Returns 0 on success, -1 on error.
 */
static int merge_collapse(mergestate_t *ms)
{
    struct s_slice *p = ms->pending;

    assert(ms);
    while (ms->n > 1) {
        signed_size_t n = ms->n - 2;
        if ((n > 0 && p[n-1].len <= p[n].len + p[n+1].len) ||
            (n > 1 && p[n-2].len <= p[n-1].len + p[n].len)) {
            if (p[n-1].len < p[n+1].len)
                --n;
            if (merge_at(ms, n) < 0)
                return -1;
        }
        else if (p[n].len <= p[n+1].len) {
            if (merge_at(ms, n) < 0)
                return -1;
        }
        else
            break;
    }
    return 0;
}

static list_t * list_t_sort_impl(list_t *self, int reverse) {

    mergestate_t ms;
    signed_size_t nremaining;
    signed_size_t minrun;
    sortslice_t lo;
    signed_size_t saved_ob_size;
    item_t *saved_ob_item;
    item_t *final_ob_item;
    list_t *result = NULL;            /* guilty until proved innocent */

    // pass the reference to the context we are in, in particular to access the Lua state.
    ms.list_tobject = self;

    /* The list_t is temporarily made empty, so that mutations performed
     * by comparison functions can't affect the slice of memory we're
     * sorting (allowing mutations during sorting is a core-dump
     * factory, since ob_item may change).
     */
    saved_ob_size = self->ob_size;
    saved_ob_item = self->ob_item;

    self->ob_size = 0;
    self->ob_item = NULL;
    
    lo  = saved_ob_item;

    merge_init(&ms, saved_ob_size);

    nremaining = saved_ob_size;
    if (nremaining < 2)
        goto succeed;

    /* Reverse sort stability achieved by initially reversing the list_t,
    applying a stable forward sort, then reversing the final result. */
    if (reverse) {
        reverse_slice(&saved_ob_item[0], &saved_ob_item[saved_ob_size]);
    }

    /* March over the array once, left to right, finding natural runs,
     * and extending short natural runs to minrun elements.
     */
    minrun = merge_compute_minrun(nremaining);
    do {
        int descending;
        signed_size_t n;

        /* Identify next run. */
        n = count_run(&ms, lo , lo  + nremaining, &descending);
        if (n < 0)
            goto fail;
        if (descending)
            reverse_sortslice_t(&lo, n);
        /* If short, extend to min(minrun, nremaining). */
        if (n < minrun) {
            const signed_size_t force = nremaining <= minrun ? nremaining : minrun;
            if (binarysort(&ms, lo, lo  + force, lo  + n) < 0)
                goto fail;
            n = force;
        }
        /* Push run onto pending-runs stack, and maybe merge. */
        assert(ms.n < MAX_MERGE_PENDING);
        ms.pending[ms.n].base = lo;
        ms.pending[ms.n].len = n;
        ++ms.n;
        if (merge_collapse(&ms) < 0)
            goto fail;
        /* Advance to find next run. */
        sortslice_t_advance(&lo, n);
        nremaining -= n;
    } while (nremaining);

    if (merge_force_collapse(&ms) < 0)
        goto fail;
    assert(ms.n == 1);
    assert(ms.pending[0].base  == saved_ob_item);
    assert(ms.pending[0].len == saved_ob_size);
    lo = ms.pending[0].base;

succeed:
    result = self;
fail:

    if (reverse && saved_ob_size > 1)
        reverse_slice(saved_ob_item, saved_ob_item + saved_ob_size);

    merge_freemem(&ms);

//keyfunc_fail:
    final_ob_item = self->ob_item;
    self->ob_size = saved_ob_size;
    self->ob_item = saved_ob_item;
    if (final_ob_item != NULL) {
        /* we cannot use _list_t_clear() for this because it does not
           guarantee that the list_t is really empty when it returns */
        /*while (--i >= 0) {
            Py_XDECREF(final_ob_item[i]);
        }*/
        free(final_ob_item);
    }

    return result;    
}

/* An adaptive, stable, natural mergesort.  See list_tsort.txt.
 * Returns Py_None on success, NULL on error.  Even in case of error, the
 * list_t will be some permutation of its input state (nothing is lost or
 * duplicated).
 */
/*[clinic input]
list_t.sort

    *
    key as keyfunc: object = None
    reverse: bool(accept={int}) = False

Sort the list_t in ascending order and return None.

The sort is in-place (i.e. the list_t itself is modified) and stable (i.e. the
order of two equal elements is maintained).

If a key function is given, apply it once to each list_t item and sort them,
ascending or descending, according to their function values.

The reverse flag can be set to sort in descending order.
[clinic start generated code]*/
static int l_sort(lua_State *L) {
	
	list_t self;

    //int nargs = lua_gettop(L);

	/* Initial checks */
	assert(lua_istable(L, -3));
	assert(lua_isboolean(L, -2));
    assert(lua_isfunction(L, -1));

    self.L = L;
    self.table_absidx = lua_absindex(L, -3);

    lua_len(L, self.table_absidx);  // push on the stack the number of elements to sort.
    int nel = lua_tointeger(L, -1); // get that number.
    lua_pop(L, 1);                  // clean the stack.

    self.ob_item = (item_t *) malloc (sizeof(item_t ) * nel);
    self.ob_size = nel;

    for (int i = 0; i < nel; i++) {
        self.ob_item[i] = i + 1;   // simply prepare the identity permutation.
    }

    int reverse = lua_toboolean(L, -2);

    time_t starttime = time(NULL);
    list_t *result = list_t_sort_impl(&self, reverse);
    time_t endtime = time(NULL);

    if(result == NULL) {
        lua_pushstring(L, "timsort error");
        lua_error(L);
    }

    lua_createtable(L, nel, 0);	// the sorted table
    lua_createtable(L, nel, 0);	// the sorting permutation

    for (int i = 0; i < nel; i++) {
        int idx = self.ob_item[i];

        lua_geti(L, self.table_absidx, idx);
        lua_seti(L, -3, i + 1);

	lua_pushinteger(L, idx);	// to also provide the sorting permutation.
        lua_seti(L, -2, i + 1);
    }

    free(self.ob_item); // finally release the all vector.

    lua_pushinteger(L, endtime - starttime);

    return 3;
}

static const struct luaL_Reg timsort [] = {
	{"sort", l_sort},
	{NULL, NULL} /* sentinel */
};
 
int luaopen_libtimsort (lua_State *L) {
	luaL_newlib(L, timsort);
	return 1;
}
