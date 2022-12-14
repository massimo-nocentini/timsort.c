
/* 
	This is a glue c file for importing delta client c functions into Lua workflow.
*/

#include <stdint.h>


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

typedef int (*comparer_t)(item_t v, item_t w, mergestate_t *ms);

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

    comparer_t cmp;

    void *userdata;
};

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
sortslice_t timsort (sortslice_t, signed_size_t, int, comparer_t, void *);
