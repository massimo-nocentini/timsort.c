

#include <limits.h>

#define PY_SSIZE_T_MAX __LONG_MAX__
#define Py_LOCAL_INLINE(type) static inline type
#define SIZEOF_SIZE_T sizeof(size_t)
#define MAX_MERGE_PENDING (SIZEOF_SIZE_T * 8)

/* The maximum number of entries in a MergeState's pending-runs stack.
 * For a list with n elements, this needs at most floor(log2(n)) + 1 entries
 * even if we didn't force runs to a minimal length.  So the number of bits
 * in a Py_ssize_t is plenty large enough for all cases.
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
// #if ((MAX_MINRUN) < 1) || ((MAX_MINRUN) & ((MAX_MINRUN) - 1))
// #error "MAX_MINRUN must be a power of 2, and >= 1"
// #endif

typedef __ssize_t ssize_t;
typedef ssize_t Py_ssize_t;

typedef struct PyObject_s
{
    void *external_object;
    int index;
} PyObject;

/* Lots of code for an adaptive, stable, natural mergesort.  There are many
 * pieces to this algorithm; read listsort.txt for overviews and details.
 */

/* A sortslice contains a pointer to an array of keys and a pointer to
 * an array of corresponding values.  In other words, keys[i]
 * corresponds with values[i].  If values == NULL, then the keys are
 * also the values.
 *
 * Several convenience routines are provided here, so that keys and
 * values are always moved in sync.
 */

typedef struct
{
    PyObject **keys;
    PyObject **values;
} sortslice;

/* One MergeState exists on the stack per invocation of mergesort.  It's just
 * a convenient way to pass state around among the helper functions.
 */
struct s_slice
{
    sortslice base;
    Py_ssize_t len; /* length of run */
    int power;      /* node "level" for powersort merge strategy */
};

typedef struct s_MergeState MergeState;

struct s_MergeState
{
    /* This controls when we get *into* galloping mode.  It's initialized
     * to MIN_GALLOP.  merge_lo and merge_hi tend to nudge it higher for
     * random data, and lower for highly structured data.
     */
    Py_ssize_t min_gallop;

    Py_ssize_t listlen;  /* len(input_list) - read only */
    PyObject **basekeys; /* base address of keys array - read only */

    /* 'a' is temp storage to help with merges.  It contains room for
     * alloced entries.
     */
    sortslice a; /* may point to temparray below */
    Py_ssize_t alloced;

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
    PyObject *temparray[MERGESTATE_TEMP_SIZE];

    /* This is the function we will use to compare two keys,
     * even when none of our special cases apply and we have to use
     * safe_object_compare. */
    int (*key_compare)(PyObject *, PyObject *, MergeState *);

    /* This function is used by unsafe_object_compare to optimize comparisons
     * when we know our list is type-homogeneous but we can't assume anything else.
     * In the pre-sort check it is set equal to Py_TYPE(key)->tp_richcompare */
    PyObject *(*key_richcompare)(PyObject *, PyObject *, int);

    /* This function is used by unsafe_tuple_compare to compare the first elements
     * of tuples. It may be set to safe_object_compare, but the idea is that hopefully
     * we can assume more, and use one of the special-case compares. */
    int (*tuple_elem_compare)(PyObject *, PyObject *, MergeState *);
};

#define Py_SIZE(ob) ob->ob_base.ob_size
#define Py_SET_SIZE(ob, size) ob->ob_base.ob_size = size
#define PyMem_FREE(p) free((p))

typedef struct
{
    PyObject ob_base;
    Py_ssize_t ob_size; /* Number of items in variable part */
} PyVarObject;

typedef struct PyListObject_s
{
    PyVarObject ob_base;
    PyObject **ob_item;
    Py_ssize_t allocated;
} PyListObject;

typedef PyObject *(*keyfunc_t)(PyObject *obj, MergeState *ms);

PyObject *list_sort_impl(PyListObject *self, keyfunc_t keyfunc, int reverse);