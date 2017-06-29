OBJ=\
bitmask.o\
chooser.o\
graph.o\
io.o\
params.o\
precode.o\
rand.o\
nanorq.o

CPPFLAGS = -D_DEFAULT_SOURCE -D_FILE_OFFSET_BITS=64 
CFLAGS   = -O2 -g -std=c99 -Wall -funroll-loops -I. -Ioblas
#LDFLAGS+= -lprofiler

all: test libnanorq.a

test: encode decode
	$(MAKE) -f example.make

encode: encode.o libnanorq.a

decode: decode.o libnanorq.a

oblas/liboblas.a:
	$(MAKE) -C oblas

.PHONY: oblas_clean
oblas_clean:
	$(MAKE) -C oblas clean

libnanorq.a: $(OBJ) oblas/liboblas.a
	$(AR) rcs $@ $^ oblas/octmat.o oblas/oblas.o

clean: oblas_clean
	$(RM) encode decode *.o *.a

indent:
	clang-format -style=LLVM -i *.c *.h

scan:
	scan-build $(MAKE) clean all
