#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#include "io.h"

struct fileioctx {
  struct ioctx io;
  FILE *fp;
};

static size_t fileio_read(struct ioctx *io, void *buf, int len) {
  struct fileioctx *fio = (struct fileioctx *)io;
  return fread(buf, 1, len, fio->fp);
}

static size_t fileio_write(struct ioctx *io, const void *buf, int len) {
  struct fileioctx *fio = (struct fileioctx *)io;
  return fwrite(buf, 1, len, fio->fp);
}

static int fileio_seek(struct ioctx *io, const int offset) {
  struct fileioctx *fio = (struct fileioctx *)io;
  return (fseek(fio->fp, offset, SEEK_SET) == 0);
}

static long fileio_tell(struct ioctx *io) {
  struct fileioctx *fio = (struct fileioctx *)io;
  return ftell(fio->fp);
}

static void fileio_destroy(struct ioctx *io) {
  struct fileioctx *fio = (struct fileioctx *)io;
  fclose(fio->fp);
  free(fio);
  return;
}

static size_t fileio_size(struct ioctx *io) {
  struct fileioctx *fio = (struct fileioctx *)io;
  long ret = 0;
  long pos = ftell(fio->fp);
  fseek(fio->fp, 0, SEEK_END);
  ret = ftell(fio->fp);
  fseek(fio->fp, pos, SEEK_SET);
  return ret;
}

struct ioctx *ioctx_from_file(const char *fn, int t) {
  struct fileioctx *ret = NULL;

  FILE *fp;

  if (t) {
    fp = fopen(fn, "r");
  } else {
    fp = fopen(fn, "w+"); // create decoder
  }

  if (!fp)
    return NULL;

  ret = calloc(1, sizeof(struct fileioctx));
  ret->fp = fp;

  ret->io.read = fileio_read;
  ret->io.write = fileio_write;
  ret->io.seek = fileio_seek;
  ret->io.size = fileio_size;
  ret->io.tell = fileio_tell;
  ret->io.destroy = fileio_destroy;
  ret->io.seekable = true;

  return (struct ioctx *)ret;
}
