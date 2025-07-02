
/*
    This is a glue c file for importing delta client c functions into Lua workflow.
*/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "timsort.h"

/* Reverse a slice of a list in place, from lo up to (exclusive) hi. */
static void
reverse_slice(timsort_object_t **lo, timsort_object_t **hi)
{
    assert(lo && hi);

    --hi;
    while (lo < hi)
    {
        timsort_object_t *t = *lo;
        *lo = *hi;
        *hi = t;
        ++lo;
        --hi;
    }
}

TIMSORT_LOCAL_INLINE(void)
sortslice_copy(timsort_sortslice_t *s1, timsort_ssize_t i, timsort_sortslice_t *s2, timsort_ssize_t j)
{
    s1->keys[i] = s2->keys[j];
}

TIMSORT_LOCAL_INLINE(void)
sortslice_copy_incr(timsort_sortslice_t *dst, timsort_sortslice_t *src)
{
    *dst->keys++ = *src->keys++;
}

TIMSORT_LOCAL_INLINE(void)
sortslice_copy_decr(timsort_sortslice_t *dst, timsort_sortslice_t *src)
{
    *dst->keys-- = *src->keys--;
}

TIMSORT_LOCAL_INLINE(void)
sortslice_memcpy(timsort_sortslice_t *s1, timsort_ssize_t i, timsort_sortslice_t *s2, timsort_ssize_t j,
                 timsort_ssize_t n)
{
    memcpy(&s1->keys[i], &s2->keys[j], sizeof(timsort_object_t *) * n);
}

TIMSORT_LOCAL_INLINE(void)
sortslice_memmove(timsort_sortslice_t *s1, timsort_ssize_t i, timsort_sortslice_t *s2, timsort_ssize_t j,
                  timsort_ssize_t n)
{
    memmove(&s1->keys[i], &s2->keys[j], sizeof(timsort_object_t *) * n);
}

TIMSORT_LOCAL_INLINE(void)
sortslice_advance(timsort_sortslice_t *slice, timsort_ssize_t n)
{
    slice->keys += n;
}

/* Comparison function: ms->key_compare, which is set at run-time in
 * listsort_impl to optimize for various special cases.
 * Returns -1 on error, 1 if x < y, 0 if x >= y.
 */

#define ISLT(X, Y) (*(ms->key_compare))(X, Y, ms)

/* Compare X to Y via "<".  Goto "fail" if the comparison raises an
   error.  Else "k" is set to true iff X<Y, and an "if (k)" block is
   started.  It makes more sense in context <wink>.  X and Y are timsort_object_t*s.
*/
#define IFLT(X, Y)            \
    if ((k = ISLT(X, Y)) < 0) \
        goto fail;            \
    if (k)

/* binarysort is the best method for sorting small arrays: it does few
   compares, but can do data movement quadratic in the number of elements.
   ss->keys is viewed as an array of n kays, a[:n]. a[:ok] is already sorted.
   Pass ok = 0 (or 1) if you don't know.
   It's sorted in-place, by a stable binary insertion sort. If ss->values
   isn't NULL, it's permuted in lockstap with ss->keys.
   On entry, must have n >= 1, and 0 <= ok <= n <= MAX_MINRUN.
   Return -1 if comparison raises an exception, else 0.
   Even in case of error, the output slice will be some permutation of
   the input (nothing is lost or duplicated).
*/
static int
binarysort(timsort_mergestate_t *ms, const timsort_sortslice_t *ss, timsort_ssize_t n, timsort_ssize_t ok)
{
    timsort_ssize_t k; /* for IFLT macro expansion */
    timsort_object_t **const a = ss->keys;

    timsort_object_t *pivot;
    timsort_ssize_t M;

    assert(0 <= ok && ok <= n && 1 <= n && n <= MAX_MINRUN);
    /* assert a[:ok] is sorted */
    if (!ok)
        ++ok;

    /* Regular insertion sort has average- and worst-case O(n**2) cost
       for both # of comparisons and number of bytes moved. But its branches
       are highly predictable, and it loves sorted input (n-1 compares and no
       data movement). This is significant in cases like sortperf.py's %sort,
       where an out-of-order element near the start of a run is moved into
       place slowly but then the remaining elements up to length minrun are
       generally at worst one slot away from their correct position (so only
       need 1 or 2 commpares to resolve). If comparisons are very fast (such
       as for a list of Python floats), the simple inner loop leaves it
       very competitive with binary insertion, despite that it does
       significantly more compares overall on random data.

       Binary insertion sort has worst, average, and best case O(n log n)
       cost for # of comparisons, but worst and average case O(n**2) cost
       for data movement. The more expensive comparisons, the more important
       the comparison advantage. But its branches are less predictable the
       more "randomish" the data, and that's so significant its worst case
       in real life is random input rather than reverse-ordered (which does
       about twice the data movement than random input does).

       Note that the number of bytes moved doesn't seem to matter. MAX_MINRUN
       of 64 is so small that the key and value pointers all fit in a corner
       of L1 cache, and moving things around in that is very fast. */
    if (ms->use_ordinary_insertion_sort)
    { // ordinary insertion sort.
        for (; ok < n; ++ok)
        {
            pivot = a[ok];

            for (M = ok - 1; M >= 0; --M)
            {
                k = ISLT(pivot, a[M]);
                if (k < 0)
                {
                    a[M + 1] = pivot;
                    goto fail;
                }
                else if (k)
                {
                    a[M + 1] = a[M];
                }
                else
                    break;
            }
            a[M + 1] = pivot;
        }
    }
    else
    { // binary insertion sort
        timsort_ssize_t L, R;
        for (; ok < n; ++ok)
        {
            /* set L to where a[ok] belongs */
            L = 0;
            R = ok;
            pivot = a[ok];
            /* Slice invariants. vacuously true at the start:
             * all a[0:L]  <= pivot
             * all a[L:R]     unknown
             * all a[R:ok]  > pivot
             */
            assert(L < R);
            do
            {
                /* don't do silly ;-) things to prevent overflow when finding
                   the midpoint; L and R are very far from filling a timsort_ssize_t */
                M = (L + R) >> 1;

                if (ms->unpredictable_branch_on_random_data)
                { // straightforward, but highly unpredictable branch on random data
                    IFLT(pivot, a[M])
                    R = M;
                    else L = M + 1;
                }
                else
                {
                    /* Try to get compiler to generate conditional move instructions
                       instead. Works fine, but leaving it disabled for now because
                       it's not yielding consistently faster sorts. Needs more
                       investigation. More computation in the inner loop adds its own
                       costs, which can be significant when compares are fast. */

                    k = ISLT(pivot, a[M]);
                    if (k < 0)
                        goto fail;
                    timsort_ssize_t Mp1 = M + 1;
                    R = k ? M : R;
                    L = k ? L : Mp1;
                }

            } while (L < R);
            assert(L == R);
            /* a[:L] holds all elements from a[:ok] <= pivot now, so pivot belongs
               at index L. Slide a[L:ok] to the right a slot to make room for it.
               Caution: using memmove is much slower under MSVC 5; we're not
               usually moving many slots. Years later: under Visual Studio 2022,
               memmove seems just slightly slower than doing it "by hand". */
            for (M = ok; M > L; --M)
                a[M] = a[M - 1];
            a[L] = pivot;
        }
    } // pick binary or regular insertion sort
    return 0;

fail:
    return -1;
}

static void
sortslice_reverse(timsort_sortslice_t *s, timsort_ssize_t n)
{
    reverse_slice(s->keys, &s->keys[n]);
}

/*
Return the length of the run beginning at slo->keys, spanning no more than
nremaining elements. The run beginning there may be ascending or descending,
but the function permutes it in place, if needed, so that it's always ascending
upon return.

Returns -1 in case of error.
*/
static timsort_ssize_t
count_run(timsort_mergestate_t *ms, timsort_sortslice_t *slo, timsort_ssize_t nremaining)
{
    timsort_ssize_t k; /* used by IFLT macro expansion */
    timsort_ssize_t n;
    timsort_object_t **const lo = slo->keys;

    /* In general, as things go on we've established that the slice starts
       with a monotone run of n elements, starting at lo. */

    /* We're n elements into the slice, and the most recent neq+1 elments are
     * all equal. This reverses them in-place, and resets neq for reuse.
     */
#define REVERSE_LAST_NEQ                    \
    if (neq)                                \
    {                                       \
        timsort_sortslice_t slice = *slo;   \
        ++neq;                              \
        sortslice_advance(&slice, n - neq); \
        sortslice_reverse(&slice, neq);     \
        neq = 0;                            \
    }

    /* Sticking to only __lt__ compares is confusing and error-prone. But in
     * this routine, almost all uses of IFLT can be captured by tiny macros
     * giving mnemonic names to the intent. Note that inline functions don't
     * work for this (IFLT expands to code including `goto fail`).
     */
#define IF_NEXT_LARGER IFLT(lo[n - 1], lo[n])
#define IF_NEXT_SMALLER IFLT(lo[n], lo[n - 1])

    assert(nremaining);
    /* try ascending run first */
    for (n = 1; n < nremaining; ++n)
    {
        IF_NEXT_SMALLER
        break;
    }
    if (n == nremaining)
        return n;
    /* lo[n] is strictly less */
    /* If n is 1 now, then the first compare established it's a descending
     * run, so fall through to the descending case. But if n > 1, there are
     * n elements in an ascending run terminated by the strictly less lo[n].
     * If the first key < lo[n-1], *somewhere* along the way the sequence
     * increased, so we're done (there is no descending run).
     * Else first key >= lo[n-1], which implies that the entire ascending run
     * consists of equal elements. In that case, this is a descending run,
     * and we reverse the all-equal prefix in-place.
     */
    if (n > 1)
    {
        IFLT(lo[0], lo[n - 1])
        return n;
        sortslice_reverse(slo, n);
    }
    ++n; /* in all cases it's been established that lo[n] has been resolved */

    /* Finish descending run. All-squal subruns are reversed in-place on the
     * fly. Their original order will be restored at the end by the whole-slice
     * reversal.
     */
    timsort_ssize_t neq = 0;
    for (; n < nremaining; ++n)
    {
        IF_NEXT_SMALLER
        {
            /* This ends the most recent run of equal elments, but still in
             * the "descending" direction.
             */
            REVERSE_LAST_NEQ
        }
        else
        {
            IF_NEXT_LARGER /* descending run is over */
                break;
            else /* not x < y and not y < x implies x == y */
                ++neq;
        }
    }
    REVERSE_LAST_NEQ
    sortslice_reverse(slo, n); /* transform to ascending run */

    /* And after reversing, it's possible this can be extended by a
     * naturally increasing suffix; e.g., [3, 2, 3, 4, 1] makes an
     * ascending run from the first 4 elements.
     */
    for (; n < nremaining; ++n)
    {
        IF_NEXT_SMALLER
        break;
    }

    return n;
fail:
    return -1;

#undef REVERSE_LAST_NEQ
#undef IF_NEXT_SMALLER
#undef IF_NEXT_LARGER
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

Returns -1 on error.  See listsort.txt for info on the method.
*/
static timsort_ssize_t
gallop_left(timsort_mergestate_t *ms, timsort_object_t *key, timsort_object_t **a, timsort_ssize_t n, timsort_ssize_t hint)
{
    timsort_ssize_t ofs;
    timsort_ssize_t lastofs;
    timsort_ssize_t k;

    assert(key && a && n > 0 && hint >= 0 && hint < n);

    a += hint;
    lastofs = 0;
    ofs = 1;
    IFLT(*a, key)
    {
        /* a[hint] < key -- gallop right, until
         * a[hint + lastofs] < key <= a[hint + ofs]
         */
        const timsort_ssize_t maxofs = n - hint; /* &a[n-1] is highest */
        while (ofs < maxofs)
        {
            IFLT(a[ofs], key)
            {
                lastofs = ofs;
                assert(ofs <= (TIMSORT_SSIZE_T_MAX - 1) / 2);
                ofs = (ofs << 1) + 1;
            }
            else /* key <= a[hint + ofs] */
                break;
        }
        if (ofs > maxofs)
            ofs = maxofs;
        /* Translate back to offsets relative to &a[0]. */
        lastofs += hint;
        ofs += hint;
    }
    else
    {
        /* key <= a[hint] -- gallop left, until
         * a[hint - ofs] < key <= a[hint - lastofs]
         */
        const timsort_ssize_t maxofs = hint + 1; /* &a[0] is lowest */
        while (ofs < maxofs)
        {
            IFLT(*(a - ofs), key)
            break;
            /* key <= a[hint - ofs] */
            lastofs = ofs;
            assert(ofs <= (TIMSORT_SSIZE_T_MAX - 1) / 2);
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
    while (lastofs < ofs)
    {
        timsort_ssize_t m = lastofs + ((ofs - lastofs) >> 1);

        IFLT(a[m], key)
        lastofs = m + 1; /* a[m] < key */
        else ofs = m;    /* key <= a[m] */
    }
    assert(lastofs == ofs); /* so a[ofs-1] < key <= a[ofs] */
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
static timsort_ssize_t
gallop_right(timsort_mergestate_t *ms, timsort_object_t *key, timsort_object_t **a, timsort_ssize_t n, timsort_ssize_t hint)
{
    timsort_ssize_t ofs;
    timsort_ssize_t lastofs;
    timsort_ssize_t k;

    assert(key && a && n > 0 && hint >= 0 && hint < n);

    a += hint;
    lastofs = 0;
    ofs = 1;
    IFLT(key, *a)
    {
        /* key < a[hint] -- gallop left, until
         * a[hint - ofs] <= key < a[hint - lastofs]
         */
        const timsort_ssize_t maxofs = hint + 1; /* &a[0] is lowest */
        while (ofs < maxofs)
        {
            IFLT(key, *(a - ofs))
            {
                lastofs = ofs;
                assert(ofs <= (TIMSORT_SSIZE_T_MAX - 1) / 2);
                ofs = (ofs << 1) + 1;
            }
            else /* a[hint - ofs] <= key */
                break;
        }
        if (ofs > maxofs)
            ofs = maxofs;
        /* Translate back to positive offsets relative to &a[0]. */
        k = lastofs;
        lastofs = hint - ofs;
        ofs = hint - k;
    }
    else
    {
        /* a[hint] <= key -- gallop right, until
         * a[hint + lastofs] <= key < a[hint + ofs]
         */
        const timsort_ssize_t maxofs = n - hint; /* &a[n-1] is highest */
        while (ofs < maxofs)
        {
            IFLT(key, a[ofs])
            break;
            /* a[hint + ofs] <= key */
            lastofs = ofs;
            assert(ofs <= (TIMSORT_SSIZE_T_MAX - 1) / 2);
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
    while (lastofs < ofs)
    {
        timsort_ssize_t m = lastofs + ((ofs - lastofs) >> 1);

        IFLT(key, a[m])
        ofs = m;              /* key < a[m] */
        else lastofs = m + 1; /* a[m] <= key */
    }
    assert(lastofs == ofs); /* so a[ofs-1] <= key < a[ofs] */
    return ofs;

fail:
    return -1;
}

/* Conceptually a timsort_mergestate_t's constructor. */
static void
merge_init(timsort_mergestate_t *ms, timsort_ssize_t list_size, timsort_sortslice_t *lo)
{
    assert(ms != NULL);
    ms->alloced = MERGESTATE_TEMP_SIZE;
    ms->a.keys = ms->temparray;
    ms->n = 0;
    ms->min_gallop = MIN_GALLOP;
    ms->listlen = list_size;
    ms->basekeys = lo->keys;
}

/* Free all the temp memory owned by the timsort_mergestate_t.  This must be called
 * when you're done with a timsort_mergestate_t, and may be called before then if
 * you want to free the temp memory early.
 */
static void
merge_freemem(timsort_mergestate_t *ms)
{
    assert(ms != NULL);
    if (ms->a.keys != ms->temparray)
    {
        free(ms->a.keys);
        ms->a.keys = NULL;
    }
}

/* Ensure enough temp memory for 'need' array slots is available.
 * Returns 0 on success and -1 if the memory can't be gotten.
 */
static int
merge_getmem(timsort_mergestate_t *ms, timsort_ssize_t need)
{

    assert(ms != NULL);
    if (need <= ms->alloced)
        return 0;

    /* Don't realloc!  That can cost cycles to copy the old data, but
     * we don't care what's in the block.
     */
    merge_freemem(ms);
    if ((size_t)need > TIMSORT_SSIZE_T_MAX / sizeof(timsort_object_t *))
    {
        // PyErr_NoMemory();
        return -1;
    }
    ms->a.keys = (timsort_object_t **)malloc(need * sizeof(timsort_object_t *));
    if (ms->a.keys != NULL)
    {
        ms->alloced = need;

        return 0;
    }
    // PyErr_NoMemory();
    return -1;
}
#define MERGE_GETMEM(MS, NEED) ((NEED) <= (MS)->alloced ? 0 : merge_getmem(MS, NEED))

/* Merge the na elements starting at ssa with the nb elements starting at
 * ssb.keys = ssa.keys + na in a stable way, in-place.  na and nb must be > 0.
 * Must also have that ssa.keys[na-1] belongs at the end of the merge, and
 * should have na <= nb.  See listsort.txt for more info.  Return 0 if
 * successful, -1 if error.
 */
static timsort_ssize_t
merge_lo(timsort_mergestate_t *ms, timsort_sortslice_t ssa, timsort_ssize_t na, timsort_sortslice_t ssb, timsort_ssize_t nb)
{
    timsort_ssize_t k;
    timsort_sortslice_t dest;
    int result = -1; /* guilty until proved innocent */
    timsort_ssize_t min_gallop;

    assert(ms && ssa.keys && ssb.keys && na > 0 && nb > 0);
    assert(ssa.keys + na == ssb.keys);
    if (MERGE_GETMEM(ms, na) < 0)
        return -1;
    sortslice_memcpy(&ms->a, 0, &ssa, 0, na);
    dest = ssa;
    ssa = ms->a;

    sortslice_copy_incr(&dest, &ssb);
    --nb;
    if (nb == 0)
        goto Succeed;
    if (na == 1)
        goto CopyB;

    min_gallop = ms->min_gallop;
    for (;;)
    {
        timsort_ssize_t acount = 0; /* # of times A won in a row */
        timsort_ssize_t bcount = 0; /* # of times B won in a row */

        /* Do the straightforward thing until (if ever) one run
         * appears to win consistently.
         */
        for (;;)
        {
            assert(na > 1 && nb > 0);
            k = ISLT(ssb.keys[0], ssa.keys[0]);
            if (k)
            {
                if (k < 0)
                    goto Fail;
                sortslice_copy_incr(&dest, &ssb);
                ++bcount;
                acount = 0;
                --nb;
                if (nb == 0)
                    goto Succeed;
                if (bcount >= min_gallop)
                    break;
            }
            else
            {
                sortslice_copy_incr(&dest, &ssa);
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
        do
        {
            assert(na > 1 && nb > 0);
            min_gallop -= min_gallop > 1;
            ms->min_gallop = min_gallop;
            k = gallop_right(ms, ssb.keys[0], ssa.keys, na, 0);
            acount = k;
            if (k)
            {
                if (k < 0)
                    goto Fail;
                sortslice_memcpy(&dest, 0, &ssa, 0, k);
                sortslice_advance(&dest, k);
                sortslice_advance(&ssa, k);
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
            sortslice_copy_incr(&dest, &ssb);
            --nb;
            if (nb == 0)
                goto Succeed;

            k = gallop_left(ms, ssa.keys[0], ssb.keys, nb, 0);
            bcount = k;
            if (k)
            {
                if (k < 0)
                    goto Fail;
                sortslice_memmove(&dest, 0, &ssb, 0, k);
                sortslice_advance(&dest, k);
                sortslice_advance(&ssb, k);
                nb -= k;
                if (nb == 0)
                    goto Succeed;
            }
            sortslice_copy_incr(&dest, &ssa);
            --na;
            if (na == 1)
                goto CopyB;
        } while (acount >= MIN_GALLOP || bcount >= MIN_GALLOP);
        ++min_gallop; /* penalize it for leaving galloping mode */
        ms->min_gallop = min_gallop;
    }
Succeed:
    result = 0;
Fail:
    if (na)
        sortslice_memcpy(&dest, 0, &ssa, 0, na);
    return result;
CopyB:
    assert(na == 1 && nb > 0);
    /* The last element of ssa belongs at the end of the merge. */
    sortslice_memmove(&dest, 0, &ssb, 0, nb);
    sortslice_copy(&dest, nb, &ssa, 0);
    return 0;
}

/* Merge the na elements starting at pa with the nb elements starting at
 * ssb.keys = ssa.keys + na in a stable way, in-place.  na and nb must be > 0.
 * Must also have that ssa.keys[na-1] belongs at the end of the merge, and
 * should have na >= nb.  See listsort.txt for more info.  Return 0 if
 * successful, -1 if error.
 */
static timsort_ssize_t
merge_hi(timsort_mergestate_t *ms, timsort_sortslice_t ssa, timsort_ssize_t na,
         timsort_sortslice_t ssb, timsort_ssize_t nb)
{
    timsort_ssize_t k;
    timsort_sortslice_t dest, basea, baseb;
    int result = -1; /* guilty until proved innocent */
    timsort_ssize_t min_gallop;

    assert(ms && ssa.keys && ssb.keys && na > 0 && nb > 0);
    assert(ssa.keys + na == ssb.keys);
    if (MERGE_GETMEM(ms, nb) < 0)
        return -1;
    dest = ssb;
    sortslice_advance(&dest, nb - 1);
    sortslice_memcpy(&ms->a, 0, &ssb, 0, nb);
    basea = ssa;
    baseb = ms->a;
    ssb.keys = ms->a.keys + nb - 1;

    sortslice_advance(&ssa, na - 1);

    sortslice_copy_decr(&dest, &ssa);
    --na;
    if (na == 0)
        goto Succeed;
    if (nb == 1)
        goto CopyA;

    min_gallop = ms->min_gallop;
    for (;;)
    {
        timsort_ssize_t acount = 0; /* # of times A won in a row */
        timsort_ssize_t bcount = 0; /* # of times B won in a row */

        /* Do the straightforward thing until (if ever) one run
         * appears to win consistently.
         */
        for (;;)
        {
            assert(na > 0 && nb > 1);
            k = ISLT(ssb.keys[0], ssa.keys[0]);
            if (k)
            {
                if (k < 0)
                    goto Fail;
                sortslice_copy_decr(&dest, &ssa);
                ++acount;
                bcount = 0;
                --na;
                if (na == 0)
                    goto Succeed;
                if (acount >= min_gallop)
                    break;
            }
            else
            {
                sortslice_copy_decr(&dest, &ssb);
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
        do
        {
            assert(na > 0 && nb > 1);
            min_gallop -= min_gallop > 1;
            ms->min_gallop = min_gallop;
            k = gallop_right(ms, ssb.keys[0], basea.keys, na, na - 1);
            if (k < 0)
                goto Fail;
            k = na - k;
            acount = k;
            if (k)
            {
                sortslice_advance(&dest, -k);
                sortslice_advance(&ssa, -k);
                sortslice_memmove(&dest, 1, &ssa, 1, k);
                na -= k;
                if (na == 0)
                    goto Succeed;
            }
            sortslice_copy_decr(&dest, &ssb);
            --nb;
            if (nb == 1)
                goto CopyA;

            k = gallop_left(ms, ssa.keys[0], baseb.keys, nb, nb - 1);
            if (k < 0)
                goto Fail;
            k = nb - k;
            bcount = k;
            if (k)
            {
                sortslice_advance(&dest, -k);
                sortslice_advance(&ssb, -k);
                sortslice_memcpy(&dest, 1, &ssb, 1, k);
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
            sortslice_copy_decr(&dest, &ssa);
            --na;
            if (na == 0)
                goto Succeed;
        } while (acount >= MIN_GALLOP || bcount >= MIN_GALLOP);
        ++min_gallop; /* penalize it for leaving galloping mode */
        ms->min_gallop = min_gallop;
    }
Succeed:
    result = 0;
Fail:
    if (nb)
        sortslice_memcpy(&dest, -(nb - 1), &baseb, 0, nb);
    return result;
CopyA:
    assert(nb == 1 && na > 0);
    /* The first element of ssb belongs at the front of the merge. */
    sortslice_memmove(&dest, 1 - na, &ssa, 1 - na, na);
    sortslice_advance(&dest, -na);
    sortslice_advance(&ssa, -na);
    sortslice_copy(&dest, 0, &ssb, 0);
    return 0;
}

/* Merge the two runs at stack indices i and i+1.
 * Returns 0 on success, -1 on error.
 */
static timsort_ssize_t
merge_at(timsort_mergestate_t *ms, timsort_ssize_t i)
{
    timsort_sortslice_t ssa, ssb;
    timsort_ssize_t na, nb;
    timsort_ssize_t k;

    assert(ms != NULL);
    assert(ms->n >= 2);
    assert(i >= 0);
    assert(i == ms->n - 2 || i == ms->n - 3);

    ssa = ms->pending[i].base;
    na = ms->pending[i].len;
    ssb = ms->pending[i + 1].base;
    nb = ms->pending[i + 1].len;
    assert(na > 0 && nb > 0);
    assert(ssa.keys + na == ssb.keys);

    /* Record the length of the combined runs; if i is the 3rd-last
     * run now, also slide over the last run (which isn't involved
     * in this merge).  The current run i+1 goes away in any case.
     */
    ms->pending[i].len = na + nb;
    if (i == ms->n - 3)
        ms->pending[i + 1] = ms->pending[i + 2];
    --ms->n;

    /* Where does b start in a?  Elements in a before that can be
     * ignored (already in place).
     */
    k = gallop_right(ms, *ssb.keys, ssa.keys, na, 0);
    if (k < 0)
        return -1;
    sortslice_advance(&ssa, k);
    na -= k;
    if (na == 0)
        return 0;

    /* Where does a end in b?  Elements in b after that can be
     * ignored (already in place).
     */
    nb = gallop_left(ms, ssa.keys[na - 1], ssb.keys, nb, nb - 1);
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

/* Two adjacent runs begin at index s1. The first run has length n1, and
 * the second run (starting at index s1+n1) has length n2. The list has total
 * length n.
 * Compute the "power" of the first run. See listsort.txt for details.
 */
static int
powerloop(timsort_ssize_t s1, timsort_ssize_t n1, timsort_ssize_t n2, timsort_ssize_t n)
{
    int result = 0;
    assert(s1 >= 0);
    assert(n1 > 0 && n2 > 0);
    assert(s1 + n1 + n2 <= n);
    /* midpoints a and b:
     * a = s1 + n1/2
     * b = s1 + n1 + n2/2 = a + (n1 + n2)/2
     *
     * Those may not be integers, though, because of the "/2". So we work with
     * 2*a and 2*b instead, which are necessarily integers. It makes no
     * difference to the outcome, since the bits in the expansion of (2*i)/n
     * are merely shifted one position from those of i/n.
     */
    timsort_ssize_t a = 2 * s1 + n1; /* 2*a */
    timsort_ssize_t b = a + n1 + n2; /* 2*b */
    /* Emulate a/n and b/n one bit a time, until bits differ. */
    for (;;)
    {
        ++result;
        if (a >= n)
        { /* both quotient bits are 1 */
            assert(b >= a);
            a -= n;
            b -= n;
        }
        else if (b >= n)
        { /* a/n bit is 0, b/n bit is 1 */
            break;
        } /* else both quotient bits are 0 */
        assert(a < b && b < n);
        a <<= 1;
        b <<= 1;
    }
    return result;
}

/* The next run has been identified, of length n2.
 * If there's already a run on the stack, apply the "powersort" merge strategy:
 * compute the topmost run's "power" (depth in a conceptual binary merge tree)
 * and merge adjacent runs on the stack with greater power. See listsort.txt
 * for more info.
 *
 * It's the caller's responsibility to push the new run on the stack when this
 * returns.
 *
 * Returns 0 on success, -1 on error.
 */
static int
found_new_run(timsort_mergestate_t *ms, timsort_ssize_t n2)
{
    assert(ms);
    if (ms->n)
    {
        assert(ms->n > 0);
        struct timsort_slice_s *p = ms->pending;
        timsort_ssize_t s1 = p[ms->n - 1].base.keys - ms->basekeys; /* start index */
        timsort_ssize_t n1 = p[ms->n - 1].len;
        int power = powerloop(s1, n1, n2, ms->listlen);
        while (ms->n > 1 && p[ms->n - 2].power > power)
        {
            if (merge_at(ms, ms->n - 2) < 0)
                return -1;
        }
        assert(ms->n < 2 || p[ms->n - 2].power < power);
        p[ms->n - 1].power = power;
    }
    return 0;
}

/* Regardless of invariants, merge all runs on the stack until only one
 * remains.  This is used at the end of the mergesort.
 *
 * Returns 0 on success, -1 on error.
 */
static int
merge_force_collapse(timsort_mergestate_t *ms)
{
    struct timsort_slice_s *p = ms->pending;

    assert(ms);
    while (ms->n > 1)
    {
        timsort_ssize_t n = ms->n - 2;
        if (n > 0 && p[n - 1].len < p[n + 1].len)
            --n;
        if (merge_at(ms, n) < 0)
            return -1;
    }
    return 0;
}

/* Compute a good value for the minimum run length; natural runs shorter
 * than this are boosted artificially via binary insertion.
 *
 * If n < MAX_MINRUN return n (it's too small to bother with fancy stuff).
 * Else if n is an exact power of 2, return MAX_MINRUN / 2.
 * Else return an int k, MAX_MINRUN / 2 <= k <= MAX_MINRUN, such that n/k is
 * close to, but strictly less than, an exact power of 2.
 *
 * See listsort.txt for more info.
 */
static timsort_ssize_t
merge_compute_minrun(timsort_ssize_t n)
{
    timsort_ssize_t r = 0; /* becomes 1 if any 1 bits are shifted off */

    assert(n >= 0);
    while (n >= MAX_MINRUN)
    {
        r |= n & 1;
        n >>= 1;
    }
    return n + r;
}

/* Here we define custom comparison functions to optimize for the cases one commonly
 * encounters in practice: homogeneous lists, often of one of the basic types. */

/* This struct holds the comparison function and helper functions
 * selected in the pre-sort check. */

/* These are the special case compare functions.
 * ms->key_compare will always point to one of these: */

/* Heterogeneous compare: default, always safe to fall back on. */
static int
safe_object_compare(timsort_object_t *v, timsort_object_t *w, timsort_mergestate_t *ms)
{
    return ms->compfunc(v, w, ms->compfunc_arg);
}

/* An adaptive, stable, natural mergesort.  See listsort.txt.
 * Returns Py_None on success, NULL on error.  Even in case of error, the
 * list will be some permutation of its input state (nothing is lost or
 * duplicated).
 */

/*
Sort the list in ascending order and return None.

The sort is in-place (i.e. the list itself is modified) and stable (i.e. the
order of two equal elements is maintained).

If a key function is given, apply it once to each list item and sort them,
ascending or descending, according to their function values.

The reverse flag can be set to sort in descending order.
*/
int list_sort_impl(timsort_list_t *self,
                   int reverse,
                   int use_ordinary_insertion_sort,
                   int unpredictable_branch_on_random_data,
                   threeways_comparefunc_t compfunc,
                   void *arg)
{
    timsort_mergestate_t ms;
    timsort_ssize_t nremaining;
    timsort_ssize_t minrun;
    timsort_sortslice_t lo;
    timsort_ssize_t saved_ob_size;
    timsort_object_t **saved_ob_item;

    int result = 1; /* guilty until proved innocent */

    assert(self != NULL);

    /* The list is temporarily made empty, so that mutations performed
     * by comparison functions can't affect the slice of memory we're
     * sorting (allowing mutations during sorting is a core-dump
     * factory, since ob_item may change).
     */
    saved_ob_size = self->ob_size;
    saved_ob_item = self->ob_item;

    self->ob_size = 0;

    lo.keys = saved_ob_item;

    /* The pre-sort check: here's where we decide which compare function to use.
     * How much optimization is safe? We test for homogeneity with respect to
     * several properties that are expensive to check at compare-time, and
     * set ms appropriately. */
    if (saved_ob_size > 1)
    {
        ms.compfunc = compfunc;
        ms.compfunc_arg = arg;
        ms.key_compare = safe_object_compare;
        ms.use_ordinary_insertion_sort = use_ordinary_insertion_sort;
        ms.unpredictable_branch_on_random_data = unpredictable_branch_on_random_data;
    }
    /* End of pre-sort check: ms is now set properly! */

    merge_init(&ms, saved_ob_size, &lo);

    nremaining = saved_ob_size;
    if (nremaining < 2)
        goto succeed;

    /* Reverse sort stability achieved by initially reversing the list,
    applying a stable forward sort, then reversing the final result. */
    if (reverse)
    {

        reverse_slice(&saved_ob_item[0], &saved_ob_item[saved_ob_size]);
    }

    /* March over the array once, left to right, finding natural runs,
     * and extending short natural runs to minrun elements.
     */
    minrun = merge_compute_minrun(nremaining);
    do
    {
        timsort_ssize_t n;

        /* Identify next run. */
        n = count_run(&ms, &lo, nremaining);
        if (n < 0)
            goto fail;
        /* If short, extend to min(minrun, nremaining). */
        if (n < minrun)
        {
            const timsort_ssize_t force = nremaining <= minrun ? nremaining : minrun;
            if (binarysort(&ms, &lo, force, n) < 0)
                goto fail;
            n = force;
        }
        /* Maybe merge pending runs. */
        assert(ms.n == 0 || ms.pending[ms.n - 1].base.keys + ms.pending[ms.n - 1].len == lo.keys);

        if (found_new_run(&ms, n) < 0)
            goto fail;
        /* Push new run on stack. */
        assert(ms.n < MAX_MERGE_PENDING);
        ms.pending[ms.n].base = lo;
        ms.pending[ms.n].len = n;
        ++ms.n;
        /* Advance to find next run. */
        sortslice_advance(&lo, n);
        nremaining -= n;
    } while (nremaining);

    if (merge_force_collapse(&ms) < 0)
        goto fail;
    assert(ms.n == 1);
    assert(ms.pending[0].base.keys == saved_ob_item);
    assert(ms.pending[0].len == saved_ob_size);
    lo = ms.pending[0].base;

succeed:
    result = 0;

fail:

    if (reverse && saved_ob_size > 1)
        reverse_slice(saved_ob_item, saved_ob_item + saved_ob_size);

    merge_freemem(&ms);

    self->ob_size = saved_ob_size;

    return result;
}
#undef IFLT
#undef ISLT
