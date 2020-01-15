#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "io.h"

struct fileioctx {
  struct ioctx io;
  FILE *fp;
};

static size_t fileio_read(struct ioctx *io, void *buf, int len) {
  struct fileioctx *_io = (struct fileioctx *)io;
  return fread(buf, 1, len, _io->fp);
}

static size_t fileio_write(struct ioctx *io, const void *buf, int len) {
  struct fileioctx *_io = (struct fileioctx *)io;
  return fwrite(buf, 1, len, _io->fp);
}

static int fileio_seek(struct ioctx *io, const int offset) {
  struct fileioctx *_io = (struct fileioctx *)io;
  return (fseek(_io->fp, offset, SEEK_SET) == 0);
}

static long fileio_tell(struct ioctx *io) {
  struct fileioctx *_io = (struct fileioctx *)io;
  return ftell(_io->fp);
}

static void fileio_destroy(struct ioctx *io) {
  struct fileioctx *_io = (struct fileioctx *)io;
  fclose(_io->fp);
  free(_io);
  return;
}

static size_t fileio_size(struct ioctx *io) {
  struct fileioctx *_io = (struct fileioctx *)io;
  long ret = 0;
  long pos = ftell(_io->fp);
  fseek(_io->fp, 0, SEEK_END);
  ret = ftell(_io->fp);
  fseek(_io->fp, pos, SEEK_SET);
  return ret;
}

struct ioctx *ioctx_from_file(const char *fn, int t) {
  struct fileioctx *_io = NULL;
  FILE *fp;

  if (t) {
    fp = fopen(fn, "r");
  } else {
    fp = fopen(fn, "w+"); // create decoder
  }

  if (!fp)
    return NULL;

  _io = calloc(1, sizeof(struct fileioctx));
  _io->fp = fp;

  _io->io.read = fileio_read;
  _io->io.write = fileio_write;
  _io->io.seek = fileio_seek;
  _io->io.size = fileio_size;
  _io->io.tell = fileio_tell;
  _io->io.destroy = fileio_destroy;
  _io->io.seekable = true;

  return (struct ioctx *)_io;
}

struct memioctx {
  struct ioctx io;
  uint8_t *ptr;
  size_t pos;
  size_t size;
};

static size_t memio_read(struct ioctx *io, void *buf, int len) {
  struct memioctx *_io = (struct memioctx *)io;
  if (_io->pos + len > _io->size)
    return 0;
  memcpy(buf, _io->ptr + _io->pos, len);
  return len;
}

static size_t memio_write(struct ioctx *io, const void *buf, int len) {
  struct memioctx *_io = (struct memioctx *)io;
  if (_io->pos + len > _io->size)
    return 0;
  memcpy(_io->ptr + _io->pos, buf, len);
  return len;
}

static int memio_seek(struct ioctx *io, const int offset) {
  struct memioctx *_io = (struct memioctx *)io;
  if (offset >= _io->size)
    return false;
  _io->pos = offset;
  return true;
}

static long memio_tell(struct ioctx *io) {
  struct memioctx *_io = (struct memioctx *)io;
  return _io->pos;
}

static void memio_destroy(struct ioctx *io) {
  struct memioctx *_io = (struct memioctx *)io;
  free(_io);
  return;
}

static size_t memio_size(struct ioctx *io) {
  struct memioctx *_io = (struct memioctx *)io;
  return _io->size;
}

struct ioctx *ioctx_from_mem(const uint8_t *ptr, size_t sz) {
  struct memioctx *_io = NULL;

  _io = calloc(1, sizeof(struct memioctx));
  _io->ptr = (uint8_t *)ptr;
  _io->pos = 0;
  _io->size = sz;

  _io->io.read = memio_read;
  _io->io.write = memio_write;
  _io->io.seek = memio_seek;
  _io->io.size = memio_size;
  _io->io.tell = memio_tell;
  _io->io.destroy = memio_destroy;
  _io->io.seekable = true;

  return (struct ioctx *)_io;
}
