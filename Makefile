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

benchmark: benchmark.o libnanorq.a
	./benchmark 1280 100 5.0
	./benchmark 1280 500 5.0
	./benchmark 1280 1000 5.0
#	./benchmark 1280 5000 5.0
#	./benchmark 1280 10000 5.0

oblas/liboblas.a:
	$(MAKE) -C oblas

.PHONY: oblas_clean
oblas_clean:
	$(MAKE) -C oblas clean

libnanorq.a: $(OBJ) oblas/liboblas.a
	$(AR) rcs $@ $^ oblas/octmat.o oblas/oblas.o oblas/sparsemat.o

clean: oblas_clean
	$(RM) encode decode *.o *.a

indent:
	clang-format -style=LLVM -i *.c *.h

scan:
	scan-build $(MAKE) clean all

profile:
	$(RM) callgrind.*
	valgrind --tool=callgrind ./encode war_and_peace.txt 1280
	gprof2dot --format=callgrind callgrind.* -z main | dot -T svg > sample.svg

gprofile:
	CPUPROFILE=gperf.prof ./encode war_and_peace.txt 1280
	pprof ./encode gperf.prof --callgrind > callgrind.gperf
	gprof2dot --format=callgrind callgrind.gperf -z main | dot -T svg > sample.svg
#	pprof ./encode gperf.prof --text

