#ifndef _HEAP_H_
#define _HEAP_H_

#define HEAP_PARENT(i) ((i - 1) / 2)
#define HEAP_LEFT(i) ((i * 2) + 1)
#define HEAP_RIGHT(i) ((i * 2) + 2)

/* Generate prototypes and inline functions */
#define HEAP_PROTOTYPE(name, type)                                             \
  void name##_HEAP_SIFTUP(type *a, size_t idx);                                \
  void name##_HEAP_SIFTDOWN(type *a, size_t idx, size_t n);                    \
  void name##_HEAP_PUSH(type *a, size_t n);                                    \
  void name##_HEAP_POP(type *a, size_t n);                                     \
  void name##_HEAPIFY(type *a, size_t n);

/* Generate methods */
#define HEAP_GENERATE(name, type, cmp)                                         \
  void name##_HEAP_SIFT_UP(type *a, size_t idx) {                              \
    type val = a[idx];                                                         \
    while (idx > 0) {                                                          \
      size_t p = HEAP_PARENT(idx);                                             \
      if (cmp(val, a[p]) > 0)                                                  \
        break;                                                                 \
      a[idx] = a[p];                                                           \
      idx = p;                                                                 \
    }                                                                          \
    a[idx] = val;                                                              \
  }                                                                            \
  void name##_HEAP_SIFT_DOWN(type *a, size_t idx, size_t n) {                  \
    size_t l = idx;                                                            \
    type val = a[idx];                                                         \
    while ((l = HEAP_LEFT(l)) < n) {                                           \
      if ((l + 1 < n) && cmp(a[l], a[l + 1]) > 0)                              \
        l++;                                                                   \
      if (cmp(a[l], val) > 0)                                                  \
        break;                                                                 \
      a[idx] = a[l];                                                           \
      idx = l;                                                                 \
    }                                                                          \
    a[idx] = val;                                                              \
  }                                                                            \
  void name##_HEAP_PUSH(type *a, size_t n) { name##_HEAP_SIFT_UP(a, n - 1); }  \
  void name##_HEAP_POP(type *a, size_t n) {                                    \
    if (n < 2)                                                                 \
      return;                                                                  \
    type tmp = a[0];                                                           \
    a[0] = a[n - 1];                                                           \
    a[n - 1] = tmp;                                                            \
    name##_HEAP_SIFT_DOWN(a, 0, n - 1);                                        \
  }                                                                            \
  void name##_HEAPIFY(type *a, size_t n) {                                     \
    for (size_t idx = (n / 2) - 1; idx != (size_t)(-1); idx--) {               \
      name##_HEAP_SIFT_DOWN(a, idx, n);                                        \
    }                                                                          \
  }

#define HEAP_SIFT_UP(name, x, y) name##_HEAP_SIFT_UP(x, y)
#define HEAP_SIFT_DOWN(name, x, y, z) name##_HEAP_SIFT_DOWN(x, y, z)
#define HEAP_PUSH(name, x, y) name##_HEAP_PUSH(x, y)
#define HEAP_POP(name, x, y) name##_HEAP_POP(x, y)
#define HEAPIFY(name, x, y) name##_HEAPIFY(x, y)

#endif /* _HEAP_H_ */
