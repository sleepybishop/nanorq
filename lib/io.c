#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "io.h"

struct fileioctx {
  struct ioctx io;
  FILE *fp;
};

static size_t fileio_read(struct ioctx *io, uint8_t *buf, size_t len) {
  struct fileioctx *_io = (struct fileioctx *)io;
  return fread(buf, 1, len, _io->fp);
}

static size_t fileio_write(struct ioctx *io, const uint8_t *buf, size_t len) {
  struct fileioctx *_io = (struct fileioctx *)io;
  return fwrite(buf, 1, len, _io->fp);
}

static bool fileio_seek(struct ioctx *io, const size_t offset) {
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
  _io->io.writable = (t == 0);

  return (struct ioctx *)_io;
}

struct memioctx {
  struct ioctx io;
  uint8_t *ptr;
  size_t pos;
  size_t size;
};

static size_t memio_read(struct ioctx *io, uint8_t *buf, size_t len) {
  struct memioctx *_io = (struct memioctx *)io;
  if (_io->pos + len > _io->size) {
    size_t diff = _io->size - _io->pos;
    memcpy(buf, _io->ptr + _io->pos, diff);
    _io->pos = _io->size;
    return diff;
  }
  memcpy(buf, _io->ptr + _io->pos, len);
  _io->pos += len;
  return len;
}

static size_t memio_write(struct ioctx *io, const uint8_t *buf, size_t len) {
  struct memioctx *_io = (struct memioctx *)io;
  if (_io->pos + len > _io->size) {
    size_t diff = _io->size - _io->pos;
    memcpy(_io->ptr + _io->pos, buf, diff);
    _io->pos = _io->size;
    return diff;
  }
  memcpy(_io->ptr + _io->pos, buf, len);
  _io->pos += len;
  return len;
}

static bool memio_seek(struct ioctx *io, const size_t offset) {
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
  _io->io.writable = true;

  return (struct ioctx *)_io;
}

struct mmapioctx {
  struct ioctx io;
  int fd;
  uint8_t *ptr;
  size_t filesize;
  size_t mapsize;
  size_t offset;
  size_t pos;
  size_t lastmap;
};

static uint8_t *mmapio_mmap(size_t mapsize, bool writable, int fd,
                            size_t offset) {
  uint8_t *ptr = NULL;

  if (writable) {
    ptr = mmap(NULL, mapsize, PROT_WRITE, MAP_SHARED, fd, offset);
  } else {
    ptr = mmap(NULL, mapsize, PROT_READ, MAP_SHARED, fd, offset);
  }

  if (ptr == MAP_FAILED) {
    perror("mmap() failed: ");
    abort();
    exit(EXIT_FAILURE);
  }
  return ptr;
}

static bool mmapio_seek(struct ioctx *io, const size_t offset) {
  struct mmapioctx *_io = (struct mmapioctx *)io;

  if (offset >= _io->offset && offset < (_io->offset + _io->mapsize)) {
    _io->pos = offset;
    if (_io->io.writable) {
      _io->filesize = (offset > _io->filesize) ? offset : _io->filesize;
    }
    return true;
  }

  if (_io->io.writable) {
    if (offset < _io->offset || offset < _io->filesize) {
      munmap(_io->ptr, _io->mapsize);
      _io->offset = (offset / _io->mapsize) * _io->mapsize;
      _io->ptr =
          mmapio_mmap(_io->mapsize, _io->io.writable, _io->fd, _io->offset);
      _io->pos = offset;
      return true;
    } else {
      munmap(_io->ptr, _io->mapsize);
      _io->offset = (offset / _io->mapsize) * _io->mapsize;
      _io->filesize = _io->offset;
      ftruncate(_io->fd, _io->filesize + _io->filesize);
      _io->ptr =
          mmapio_mmap(_io->mapsize, _io->io.writable, _io->fd, _io->offset);
      _io->pos = offset;
      return true;
    }
  }

  if (!_io->io.writable && offset >= _io->filesize)
    return false;

  if (!_io->io.writable) {
    munmap(_io->ptr, _io->lastmap);
    _io->offset = (offset / _io->mapsize) * _io->mapsize;
    size_t tmp = _io->mapsize;
    if (_io->offset + tmp > _io->filesize)
      tmp = _io->filesize - _io->offset;
    _io->ptr =
        mmapio_mmap(_io->mapsize, _io->io.writable, _io->fd, _io->offset);
    _io->pos = offset;
    _io->lastmap = tmp;
    return true;
  }

  return false;
}

static size_t mmapio_read(struct ioctx *io, uint8_t *buf, size_t len) {
  struct mmapioctx *_io = (struct mmapioctx *)io;

  size_t at = _io->pos % _io->mapsize;
  if (_io->pos + len > _io->filesize) {
    size_t diff = _io->filesize - _io->pos;
    memcpy(buf, _io->ptr + at, diff);
    _io->pos = _io->filesize;
    return diff;
  }

  if (at + len < _io->mapsize) {
    memcpy(buf, _io->ptr + at, len);
    _io->pos += len;
    return len;
  } else {
    size_t read = 0, diff = 0;
    while (read < len || _io->pos == _io->filesize) {
      at = _io->pos % _io->mapsize;
      diff = (_io->lastmap - at);
      memcpy(buf, _io->ptr + at, diff);
      read += diff;
      _io->pos += diff;
      if (_io->lastmap != _io->mapsize)
        break;
      if (mmapio_seek(io, _io->offset + _io->mapsize)) {
        int left = len - diff;
        if (_io->pos + left > _io->filesize) {
          left = _io->filesize - _io->pos;
        }
        memcpy(buf + read, _io->ptr, left);
        _io->pos += left;
        read += left;
      }
    }
    return read;
  }

  return 0;
}

static size_t mmapio_write(struct ioctx *io, const uint8_t *buf, size_t len) {
  struct mmapioctx *_io = (struct mmapioctx *)io;
  size_t written = 0, at = 0;

  if (_io->pos + len < _io->offset + _io->mapsize) {
    at = _io->pos % _io->mapsize;
    memcpy(_io->ptr + at, buf + written, len);
    _io->pos += len;
    _io->filesize = (_io->pos > _io->filesize) ? _io->pos : _io->filesize;
    return len;
  }

  while (written < len) {
    at = _io->pos % _io->mapsize;
    size_t diff = (at + len < _io->mapsize) ? len : (_io->mapsize - at);
    memcpy(_io->ptr + at, buf + written, diff);
    written += diff;
    _io->pos += diff;
    _io->filesize = (_io->pos > _io->filesize) ? _io->pos : _io->filesize;
    // seek if needed
    if (mmapio_seek(io, _io->offset + _io->mapsize)) {
      size_t left = len - diff;
      at = _io->pos % _io->mapsize;
      memcpy(_io->ptr + at, buf + written, left);
      written += left;
      _io->pos += left;
      _io->filesize = (_io->pos > _io->filesize) ? _io->pos : _io->filesize;
    }
  }

  return written;
}

static long mmapio_tell(struct ioctx *io) {
  struct mmapioctx *_io = (struct mmapioctx *)io;
  return _io->pos;
}

static void mmapio_destroy(struct ioctx *io) {
  struct mmapioctx *_io = (struct mmapioctx *)io;
  munmap(_io->ptr, _io->lastmap);
  if (_io->io.writable) {
    ftruncate(_io->fd, _io->filesize);
  }
  close(_io->fd);
  free(_io);
  return;
}

static size_t mmapio_size(struct ioctx *io) {
  struct mmapioctx *_io = (struct mmapioctx *)io;
  long ret = 0;
  long pos = lseek(_io->fd, 0, SEEK_CUR);
  lseek(_io->fd, 0, SEEK_END);
  ret = lseek(_io->fd, 0, SEEK_CUR);
  lseek(_io->fd, pos, SEEK_SET);
  return ret;
}

struct ioctx *ioctx_mmap_file(const char *fn, int t) {
  struct mmapioctx *_io = NULL;
  int fd;
  uint8_t *ptr = NULL;
  size_t filesize = 0;
  size_t offset = 0;
  size_t pagesize = sysconf(_SC_PAGESIZE);
  size_t mapsize = (65536 / pagesize) * pagesize;

  if (t) {
    fd = open(fn, O_RDONLY);
  } else {
    fd = open(fn, O_RDWR | O_CREAT | O_TRUNC, 0666); // create decoder
  }

  if (fd == -1) {
    return NULL;
  }

  if (t) {
    struct stat sb;
    fstat(fd, &sb);
    filesize = sb.st_size;
    if (filesize < mapsize)
      mapsize = filesize;
    ptr = mmapio_mmap(mapsize, false, fd, offset);
  } else {
    ftruncate(fd, mapsize);
    ptr = mmapio_mmap(mapsize, true, fd, offset);
  }

  _io = calloc(1, sizeof(struct mmapioctx));
  _io->fd = fd;
  _io->ptr = ptr;
  _io->filesize = filesize;
  _io->mapsize = mapsize;
  _io->lastmap = mapsize;
  _io->offset = offset;
  _io->pos = offset;

  _io->io.read = mmapio_read;
  _io->io.write = mmapio_write;
  _io->io.seek = mmapio_seek;
  _io->io.size = mmapio_size;
  _io->io.tell = mmapio_tell;
  _io->io.destroy = mmapio_destroy;
  _io->io.seekable = true;
  _io->io.writable = (t == 0);

  return (struct ioctx *)_io;
}
