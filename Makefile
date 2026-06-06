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

TEST_UTILS=\
t/00util/matgen\
t/00util/repgen\
t/00util/hdpcgen\
t/00util/precond\
t/00util/ult\
t/00util/schedgen

EXAMPLES=\
examples/encode\
examples/decode

CPPFLAGS = -D_DEFAULT_SOURCE -D_FILE_OFFSET_BITS=64 
CFLAGS   = -O3 -g -std=c11 -Wall -I. -Iinclude -Ideps/
CFLAGS  += -march=native -funroll-loops -ftree-vectorize -fno-inline -fstack-protector-all -Wno-unused -Wno-sequence-point

all: test libnanorq.a $(EXAMPLES)

test: encode decode
	$(MAKE) -f example.make

encode: encode.o libnanorq.a 

decode: decode.o libnanorq.a

benchmark: benchmark.o $(OBJ)

benchmark_core: benchmark_core.o $(OBJ)

t/00util/matgen: t/00util/matgen.o $(OBJ)
t/00util/repgen: t/00util/repgen.o $(OBJ)
t/00util/hdpcgen: t/00util/hdpcgen.o $(OBJ)
t/00util/precond: t/00util/precond.o $(OBJ)
t/00util/ult: t/00util/ult.o $(OBJ)
t/00util/schedgen: t/00util/schedgen.o $(OBJ)

examples/encode: examples/encode.o $(OBJ)
examples/decode: examples/decode.o $(OBJ)

check: CPPFLAGS=
check: clean $(TEST_UTILS) $(EXAMPLES)
	prove -I. -v t/*.t

check-nolibc: clean
	$(MAKE) libnanorq.a CPPFLAGS="$(CPPFLAGS) -DNANORQ_NO_LIBC"
	$(MAKE) $(TEST_UTILS)
	prove -I. -v t/10pcmat.t t/20repmat.t t/30precond.t t/35ult.t t/40hdpc.t t/50schedules.t

libnanorq.a: $(OBJ)
	$(AR) rcs $@ $(OBJ)

clean:
	$(RM) encode decode lib/*.o deps/obl/*.o *.o *.a *.gcda *.gcno *.gcov callgrind.* *.gperf *.prof *.heap perf.data perf.data.old benchmark benchmark_core $(TEST_UTILS) $(EXAMPLES)
	find . -name '*.[a,o]' | xargs $(RM)

indent:
	clang-format -style=LLVM -i lib/*.c include/*.h examples/*.c t/00util/*.c benchmark.c benchmark_core.c

scan:
	scan-build $(MAKE) clean benchmark

gcov: CFLAGS += -O0 -fprofile-arcs -ftest-coverage
gcov: LDLIBS = -lgcov --coverage
gcov: clean benchmark
	./benchmark 1280 1000 5.0

perf: clean benchmark
	perf record -g ./benchmark 1280 50000 5.0
	pprof -svg ./benchmark perf.data > perf.svg

ubsan: CC=clang
ubsan: CFLAGS += -fsanitize=address,undefined,implicit-conversion,integer
ubsan: LDLIBS += -lubsan
ubsan: clean benchmark
	./benchmark 1280 50000 0

bench: benchmark
	@echo "K       encode   precalc  decode  decode-oh5"
	@./benchmark 1280  100 5.0
	@./benchmark 1280  500 5.0
	@./benchmark 1280 1000 5.0
	@./benchmark 1280 5000 5.0
	@./benchmark 1280 10000 5.0
	@./benchmark 1280 50000 5.0

bench-core: benchmark_core
	@echo "K       encode   precalc  decode  decode-oh5"
	@./benchmark_core 1280  100 5.0
	@./benchmark_core 1280  500 5.0
	@./benchmark_core 1280 1000 5.0
	@./benchmark_core 1280 5000 5.0
	@./benchmark_core 1280 10000 5.0
	@./benchmark_core 1280 50000 5.0

check-embedded:
	$(MAKE) clean
	$(MAKE) libnanorq.a CPPFLAGS="$(CPPFLAGS) -DNANORQ_NO_LIBC"
	@echo "--- Undefined symbols in libnanorq.a ---"
	@nm -u libnanorq.a | grep -E '\b(malloc|calloc|realloc|free|posix_memalign|__assert_fail)\b' && \
		(echo "FAIL: libc symbols found in embedded build" && exit 1) || \
		echo "PASS: no libc allocator/assert symbols found"

valgrind: CPPFLAGS=-Wall -Iinclude -Ideps/ -fPIC
valgrind: CFLAGS = -O0 -g -std=c11
valgrind: clean $(TEST_UTILS) $(EXAMPLES)
	valgrind --error-exitcode=2 ./t/00util/hdpcgen  500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/matgen   500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/precond  500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/repgen   500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/ult      500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/schedgen 500 > /dev/null
	valgrind --error-exitcode=2 ./examples/encode   500 64 10 t/assets/sample_data/raw > /dev/null

check-interop:
	./t/interop_harness.sh
