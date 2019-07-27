#ifndef NANORQ_IOCTX_H
#define NANORQ_IOCTX_H

#include <stdbool.h>
#include <stdint.h>

struct ioctx {
  size_t (*read)(struct ioctx *, void *, int);
  size_t (*write)(struct ioctx *, const void *, int);
  int (*seek)(struct ioctx *, const int);
  size_t (*size)(struct ioctx *);
  long (*tell)(struct ioctx *);
  void (*destroy)(struct ioctx *);
  bool seekable;
};

struct ioctx *ioctx_from_file(const char *fn, int t);
struct ioctx *ioctx_from_mem(const uint8_t *ptr, size_t t);

#endif
