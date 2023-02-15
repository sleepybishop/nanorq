OBJ=\
lib/bitmask.o\
lib/io.o\
lib/params.o\
lib/precode.o\
lib/rand.o\
lib/sched.o\
lib/spmat.o\
lib/tuple.o\
lib/wrkmat.o\
lib/nanorq.o

CPPFLAGS = -D_DEFAULT_SOURCE -D_FILE_OFFSET_BITS=64 
CFLAGS   = -O3 -g -std=c99 -Wall -I. -Iinclude -Ideps/oblas
CFLAGS  += -march=native -funroll-loops -ftree-vectorize -fno-inline -fstack-protector-all

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

deps/oblas/liboblas.a:
	$(MAKE) -C deps/oblas CPPFLAGS+="-DOBLAS_AVX -DOCTMAT_ALIGN=32"

.PHONY: oblas_clean
oblas_clean:
	$(MAKE) -C deps/oblas clean

libnanorq.a: $(OBJ) deps/oblas/liboblas.a
	$(AR) rcs $@ $(OBJ) deps/oblas/*.o

clean: oblas_clean
	$(RM) encode decode lib/*.o *.o *.a *.gcda *.gcno *.gcov callgrind.* *.gperf *.prof *.heap perf.data perf.data.old

indent:
	clang-format -style=LLVM -i *.c *.h

scan:
	scan-build $(MAKE) clean benchmark

profile:
	$(RM) callgrind.*
	valgrind --tool=callgrind ./benchmark 1280 1000 5.0
	gprof2dot --format=callgrind callgrind.out.* -z main | dot -T svg > callgrind.svg

gcov: CFLAGS += -O0 -fprofile-arcs -ftest-coverage
gcov: LDLIBS = -lgcov --coverage
gcov: clean benchmark
	./benchmark 1280 1000 5.0

gperf: LDLIBS = -lprofiler -ltcmalloc
gperf: clean benchmark
	CPUPROFILE_FREQUENCY=100000 CPUPROFILE=gperf.prof ./benchmark 1280 50000 5.0
	pprof ./benchmark gperf.prof --callgrind > callgrind.gperf
	gprof2dot --format=callgrind callgrind.gperf -z main | dot -T svg > gperf.svg
#	pprof ./encode gperf.prof --text

gheap: LDLIBS = -lprofiler -ltcmalloc
gheap: clean benchmark
	$(RM) gmem.prof.*
	HEAPPROFILE=gmem.prof HEAP_PROFILE_INUSE_INTERVAL=1024000 ./benchmark 1280 50000 5.0
	pprof ./benchmark gmem.prof.*.heap --callgrind > memgrind.gperf
	gprof2dot --format=callgrind memgrind.gperf | dot -T svg > gheap.svg

perf:
	perf record -g ./benchmark 1280 5000 5.0
	perf script | gprof2dot --format=perf | dot -T svg > perf.svg
	#perf report

ubsan: CC=clang
ubsan: CFLAGS += -fsanitize=address,undefined,implicit-conversion,integer
ubsan: LDLIBS += -lubsan
ubsan: clean benchmark
	./benchmark 1280 50000 0

