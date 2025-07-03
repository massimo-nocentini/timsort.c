
// #include <stdio.h>
// #include <limits.h>
#include <sys/types.h>

#define TIMSORT_SSIZE_T_MAX __LONG_MAX__
#define TIMSORT_LOCAL_INLINE(type) static inline type
#define SIZEOF_SIZE_T sizeof(size_t)
#define MAX_MERGE_PENDING (SIZEOF_SIZE_T * 8)

/* The maximum number of entries in a timsort_mergestate_t's pending-runs stack.
 * For a list with n elements, this needs at most floor(log2(n)) + 1 entries
 * even if we didn't force runs to a minimal length.  So the number of bits
 * in a timsort_ssize_t is plenty large enough for all cases.
 */
#define MAX_MERGE_PENDING (SIZEOF_SIZE_T * 8)

/* When we get into galloping mode, we stay there until both runs win less
 * often than MIN_GALLOP consecutive times.  See listsort.txt for more info.
 */
#define MIN_GALLOP 7

/* Avoid malloc for small temp arrays. */
#define MERGESTATE_TEMP_SIZE 256

/* The largest value of minrun. This must be a power of 2, and >= 1, so that
 * the compute_minrun() algorithm guarantees to return a result no larger than
 * this,
 */
#define MAX_MINRUN 64

typedef ssize_t timsort_ssize_t;

typedef void *timsort_object_t;

/* Lots of code for an adaptive, stable, natural mergesort.  There are many
 * pieces to this algorithm; read listsort.txt for overviews and details.
 */

/* A timsort_sortslice_t contains a pointer to an array of keys and a pointer to
 * an array of corresponding values.  In other words, keys[i]
 * corresponds with values[i].  If values == NULL, then the keys are
 * also the values.
 *
 * Several convenience routines are provided here, so that keys and
 * values are always moved in sync.
 */

typedef struct
{
    timsort_object_t **keys;
} timsort_sortslice_t;

/* One timsort_mergestate_t exists on the stack per invocation of mergesort.  It's just
 * a convenient way to pass state around among the helper functions.
 */
struct timsort_slice_s
{
    timsort_sortslice_t base;
    timsort_ssize_t len; /* length of run */
    int power;           /* node "level" for powersort merge strategy */
};

typedef int (*threeways_comparefunc_t)(timsort_object_t *a, timsort_object_t *b, void *arg);

typedef struct timsort_mergestate_s timsort_mergestate_t;

struct timsort_mergestate_s
{
    /* This controls when we get *into* galloping mode.  It's initialized
     * to MIN_GALLOP.  merge_lo and merge_hi tend to nudge it higher for
     * random data, and lower for highly structured data.
     */
    timsort_ssize_t min_gallop;

    timsort_ssize_t listlen;     /* len(input_list) - read only */
    timsort_object_t **basekeys; /* base address of keys array - read only */

    /* 'a' is temp storage to help with merges.  It contains room for
     * alloced entries.
     */
    timsort_sortslice_t a; /* may point to temparray below */
    timsort_ssize_t alloced;

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
    struct timsort_slice_s pending[MAX_MERGE_PENDING];

    /* 'a' points to this when possible, rather than muck with malloc. */
    timsort_object_t *temparray[MERGESTATE_TEMP_SIZE];

    /* This is the function we will use to compare two keys,
     * even when none of our special cases apply and we have to use
     * safe_object_compare. */
    int (*key_compare)(timsort_object_t *, timsort_object_t *, timsort_mergestate_t *);

    threeways_comparefunc_t compfunc;
    void *compfunc_arg;

    int use_ordinary_insertion_sort;         /* Use ordinary insertion sort? */
    int unpredictable_branch_on_random_data; /* Use straightforward pivot selection? */
};

typedef struct
{
    timsort_ssize_t ob_size; /* Number of items in variable part */
    timsort_object_t **ob_item;
} timsort_list_t;

int list_sort_impl(timsort_list_t *self,
                   int reverse,
                   int use_ordinary_insertion_sort,
                   int unpredictable_branch_on_random_data,
                   threeways_comparefunc_t compfunc,
                   void *arg);