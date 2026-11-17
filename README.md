# timsort.c

A standalone C port of the sort of CPython's lists (`listsort` of `Objects/listobject.c`, Tim Peters' *timsort* with the *powersort* merge policy of CPython 3.13 and 3.14). It is adaptive, stable and natural:

- it finds the runs already in the data, ascending or strictly descending, and merges them, so sorted, reversed and nearly sorted inputs cost about one comparison per element;
- equal elements keep their relative order, also in descending sorts;
- it makes as few comparisons as a merge sort on random data, which matters when comparing is expensive (a callback into another language, for example).

The algorithm is CPython's, function for function: only the key-function arrays were removed, the memory functions became `malloc`/`free`, and two compile-time switches of CPython's binary insertion sort became runtime flags. It is the C engine of the TimSort bindings of Pharo: the TimSortPlugin and the library libTimSortFFI of [massimo-nocentini/pharo-vm](https://github.com/massimo-nocentini/pharo-vm), used by the Smalltalk packages of [massimo-nocentini/timsort.st](https://github.com/massimo-nocentini/timsort.st).

## Files

- `src/timsort.h`: the types and the one entry point, `list_sort_impl`.
- `src/timsort.c`: the sort; it needs only `malloc`, `free`, `memcpy` and `memmove`, and has no global state.
- `src/Makefile`: shared libraries for Linux, MinGW and macOS, and a wasm32 build with wasi-libc.

## API

```c
#include "timsort.h"

typedef void *timsort_object_t;
typedef int (*threeways_comparefunc_t)(timsort_object_t *a, timsort_object_t *b, void *arg);
typedef struct { timsort_ssize_t ob_size; timsort_object_t **ob_item; } timsort_list_t;

int list_sort_impl(timsort_list_t *self,
                   int reverse,
                   int use_ordinary_insertion_sort,
                   int unpredictable_branch_on_random_data,
                   threeways_comparefunc_t compfunc,
                   void *arg);
```

- **What is sorted.** `self->ob_item` is an array of `self->ob_size` pointer-sized words (8 bytes on 64-bit targets), permuted in place. Any value of that width can be stored directly: pointers to records, `int64_t` values, the bit patterns of `double`s, indices into another array. Narrower elements (`float`, `int32_t`) must be widened to 8-byte words, or sorted through an array of indices.
- **The comparison.** `compfunc(a, b, arg)` gets two of the stored words themselves (not the addresses of their slots) and the `arg` given to `list_sort_impl`. Despite its type name it is a **strict less-than predicate**: answer `1` when `a` must come before `b`, `0` otherwise, and a negative number to abort the sort.
  - A non-strict predicate (`a <= b`) still sorts but loses stability.
  - A three-way compare (`-1`/`0`/`1`) does not work: its `-1` aborts the sort.
  - To use a "may precede" predicate `before(a, b)` such as `a <= b`, pass `lt(a, b) = !before(b, a)`.
- **The result.** `0` on success, `1` on any failure: the comparison aborted, or `malloc` failed. Only the caller can tell the two apart (record the reason in `arg`). After a failure the array is still a permutation of its input, in no particular order.
- **`reverse`** (non-zero): a stable descending sort. Equal elements keep their original order; this is not an ascending sort reversed afterwards.
- **`use_ordinary_insertion_sort`** and **`unpredictable_branch_on_random_data`** choose the variant of the insertion sort that extends short runs. CPython's default is `0, 1` (binary insertion, plain branch): use it unless measurements say otherwise, and always when comparisons are expensive, since ordinary insertion makes more of them. The result is the same in every mode.
- **Memory.** One merge buffer of at most `ob_size / 2` words; up to 256 words it lives on the stack, beyond that it is `malloc`ed, grown by freeing and allocating again, and freed before returning. The state of a sort is on the C stack (about 3.7 KB), so sorts on different arrays may run in parallel or nest inside a comparison.
- **Arrays of fewer than 2 elements** return `0` at once, without calling `compfunc`.

## Example

```c
#include <stdint.h>
#include "timsort.h"

static int int64_lt(timsort_object_t *a, timsort_object_t *b, void *arg)
{
    (void)arg;
    return (int64_t)(intptr_t)a < (int64_t)(intptr_t)b;
}

int sort_int64(int64_t *values, timsort_ssize_t n)
{
    timsort_list_t list = { n, (timsort_object_t **)values };
    return list_sort_impl(&list, 0, 0, 1, int64_lt, NULL);
}
```

## Building

```sh
cd src
make linux      # libtimsort.so      (make install: /usr/local/{include,lib})
make macos      # libtimsort.dylib   (make install-macos)
make mingw      # libtimsort.dll     (make install-mingw)
make wasm-32    # web/timsort.wasm, with wasi-libc (make wasi-libc installs it)
```

Compile `timsort.c` as C, not C++ (it jumps over initializations, which C++ refuses); the header can be included from C++. The asserts are cheap invariants; `-DNDEBUG` (as in the Makefile) removes them.

## Speed

Measured on an Intel Xeon Gold 6238R with clang 18 `-O2`, 8-byte entries compared through a function pointer, nanoseconds per element, 10^6 elements:

| input | timsort | glibc `qsort` | plain merge sort |
|---|---:|---:|---:|
| random `int64_t` | 136 | 103 | 72 |
| sorted | 2.4 | 27 | 39 |
| reversed | 2.8 | 31 | 44 |
| nearly sorted (1 % random) | 11 | 40 | 60 |
| few distinct values | 49 | 52 | 68 |

On random data with cheap comparisons a plain merge sort is faster. timsort wins on any input with structure, and makes about as few comparisons as possible everywhere: 1 per element on sorted or reversed input, against about 10 for `qsort`.

## Changes from the first port

Found while binding it to Pharo, with a harness of 44,000 sorts that passes natively, under AddressSanitizer/UBSan with asserts on, and in wasm64 (Emscripten `-m64`):

1. **Asserts on zero values.** `gallop_left` and `gallop_right` asserted that the key is not `NULL`. In CPython keys are never `NULL`, but here a stored value can be `0` (the integer 0, the double `+0.0`, the index 0), so a build without `-DNDEBUG` aborted. The assert no longer tests the key.
2. **64-bit Windows.** `TIMSORT_SSIZE_T_MAX` was `__LONG_MAX__`, which is 2^31-1 where `long` is 32-bit (MinGW, LLP64), while `ssize_t` is 64-bit: a merge needing more than 268,435,455 buffer words failed. `timsort_ssize_t` is now `ptrdiff_t` and its maximum `PTRDIFF_MAX` (from `<stddef.h>` and `<stdint.h>`, so the POSIX `<sys/types.h>` is no longer needed).
3. **The header** has an include guard and `extern "C"`, defines `MAX_MERGE_PENDING` once, and no longer defines `SIZEOF_SIZE_T`, a name that some `config.h` files define too.
4. **Speed.** `ISLT` calls `compfunc` directly, instead of through `safe_object_compare`: about 10 % faster.

## Licence

This code is derived from CPython's `Objects/listobject.c`, Copyright (c) 2001 Python Software Foundation, distributed under the Python Software Foundation License Version 2 (PSF-2.0, https://docs.python.org/3/license.html). This repository has no licence file of its own yet; until it has one, the PSF-2.0 terms of the original code apply to it.
