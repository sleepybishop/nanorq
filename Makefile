OBJ=\
lib/chooser.o\
lib/io.o\
lib/nanorq.o\
lib/nanorq_core.o\
lib/ops.o\
lib/params.o\
lib/precode.o\
lib/rand.o\
lib/tuple.o\
lib/uvec.o\
deps/obl/oblas_lite.o

CPPFLAGS = -D_DEFAULT_SOURCE -D_FILE_OFFSET_BITS=64 
CFLAGS   = -O3 -g -std=c11 -Wall -I. -Iinclude -Ideps/
CFLAGS  += -march=native -funroll-loops -ftree-vectorize -fno-inline -fstack-protector-all -Wno-unused -Wno-sequence-point

all: test libnanorq.a

test: encode decode
	$(MAKE) -f example.make

encode: encode.o libnanorq.a 

decode: decode.o libnanorq.a

benchmark: benchmark.c libnanorq.a
	$(CC) $(CFLAGS) benchmark.c -o $@ libnanorq.a $(LDLIBS)

bench: graph.dat
	cat graph.dat

bench.md: graph.dat
	cat graph.dat | awk -f graph.awk 

graph.dat: benchmark
	echo "K       encode   precalc  decode  decode-oh5" > graph.dat
	./benchmark 1280   100 5.0 >> graph.dat 
	./benchmark 1280   500 5.0 >> graph.dat
	./benchmark 1280  1000 5.0 >> graph.dat
	./benchmark 1280  5000 5.0 >> graph.dat	
	./benchmark 1280 10000 5.0 >> graph.dat
	./benchmark 1280 50000 5.0 >> graph.dat

graph.png: graph.dat graph.gnuplot
	gnuplot -e "argtitle='Throughput (packet size=1280) `lscpu|grep -i 'model name'|cut -f2 -d:|xargs`'" graph.gnuplot 

libnanorq.a: $(OBJ)
	$(AR) rcs $@ $(OBJ)

clean:
	$(RM) encode decode lib/*.o deps/obl/*.o *.o *.a *.gcda *.gcno *.gcov callgrind.* *.gperf *.prof *.heap perf.data perf.data.old benchmark

indent:
	clang-format -style=LLVM -i lib/*.c include/*.h

scan:
	scan-build $(MAKE) clean benchmark

gcov: CFLAGS += -O0 -fprofile-arcs -ftest-coverage
gcov: LDLIBS = -lgcov --coverage
gcov: clean benchmark
	./benchmark 1280 1000 5.0

perf: clean benchmark
	perf record -g ./benchmark 1280 50000 5.0
	pprof -svg ./benchmark perf.data > perf.svg
	#pprof ./benchmark perf.data --text

ubsan: CC=clang
ubsan: CFLAGS += -fsanitize=address,undefined,implicit-conversion,integer
ubsan: LDLIBS += -lubsan
ubsan: clean benchmark
	./benchmark 1280 50000 0
