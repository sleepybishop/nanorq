CORE_OBJ=\
lib/chooser.o\
lib/nanorq_core.o\
lib/ops.o\
lib/params.o\
lib/partition.o\
lib/sopi.o\
lib/precode.o\
lib/rand.o\
lib/tuple.o\
lib/uvec.o\
deps/obl/oblas_common.o\
deps/obl/oblas_lite.o

OBJ=\
lib/io.o\
lib/nanorq.o\
$(CORE_OBJ)


TEST_UTILS=\
t/00util/matgen\
t/00util/repgen\
t/00util/hdpcgen\
t/00util/precond\
t/00util/ult\
t/00util/schedgen\
t/00util/test_utils\
t/00util/api_regress

EXAMPLES=\
examples/encode\
examples/decode\
examples/carousel\
examples/tsnc_multipath\
examples/tsnc_sync\
examples/blockchain_gossip

CPPFLAGS = -D_DEFAULT_SOURCE -D_FILE_OFFSET_BITS=64 
CFLAGS   = -O3 -g -std=c11 -Wall -I. -Iinclude -Ideps/
CFLAGS  += -march=native -funroll-loops -ftree-vectorize -fstack-protector-strong -Wno-unused -Wno-sequence-point

all: test libnanorq.a libnanorq_core.a $(EXAMPLES)

test: encode decode
	$(MAKE) -f example.make

encode: encode.o libnanorq.a 

decode: decode.o libnanorq.a

benchmark: benchmark.o $(OBJ)

benchmark_core: benchmark_core.o $(CORE_OBJ)

t/00util/matgen: t/00util/matgen.o $(CORE_OBJ)
t/00util/repgen: t/00util/repgen.o $(CORE_OBJ)
t/00util/hdpcgen: t/00util/hdpcgen.o $(CORE_OBJ)
t/00util/precond: t/00util/precond.o $(CORE_OBJ)
t/00util/ult: t/00util/ult.o $(CORE_OBJ)
t/00util/schedgen: t/00util/schedgen.o $(CORE_OBJ)
t/00util/test_utils: t/00util/test_utils.o $(CORE_OBJ)
t/00util/api_regress: t/00util/api_regress.o $(OBJ)

examples/encode: examples/encode.o $(OBJ)
examples/decode: examples/decode.o $(OBJ)
examples/carousel: examples/carousel.o $(OBJ)
examples/tsnc_multipath: examples/tsnc_multipath.o $(OBJ)
examples/tsnc_sync: examples/tsnc_sync.o $(OBJ)
examples/blockchain_gossip: examples/blockchain_gossip.o $(OBJ)

check: clean $(TEST_UTILS) $(EXAMPLES)
	prove -I. -v t/*.t

check-nolibc: clean
	$(MAKE) libnanorq_core.a CPPFLAGS="$(CPPFLAGS) -DNANORQ_NO_LIBC"
	$(MAKE) $(TEST_UTILS)
	prove -I. -v t/10pcmat.t t/20repmat.t t/30precond.t t/35ult.t t/40hdpc.t t/50schedules.t

libnanorq.a: $(OBJ)
	$(AR) rcs $@ $(OBJ)

libnanorq_core.a: $(CORE_OBJ)
	$(AR) rcs $@ $(CORE_OBJ)

clean:
	$(RM) encode decode lib/*.o deps/obl/*.o *.o *.a *.gcda *.gcno *.gcov callgrind.* *.gperf *.prof *.heap perf.data perf.data.old benchmark benchmark_core $(TEST_UTILS) $(EXAMPLES)
	find . -name '*.[a,o]' | xargs $(RM)

indent:
	clang-format -style=LLVM -i lib/*.c include/*.h examples/*.c t/00util/*.c benchmark.c benchmark_core.c

scan:
	scan-build --status-bugs $(MAKE) clean benchmark

gcov: CFLAGS += -O0 -fprofile-arcs -ftest-coverage
gcov: LDLIBS = -lgcov --coverage
gcov: clean benchmark
	./benchmark 1280 1000 5.0

perf: clean benchmark
	perf record -g ./benchmark 1280 50000 5.0
	pprof -svg ./benchmark perf.data > perf.svg

gperf: LDLIBS += -lprofiler
gperf: clean benchmark
	CPUPROFILE=benchmark.prof ./benchmark 1280 50000 5.0
	pprof -svg ./benchmark benchmark.prof > gperf.svg

ubsan: CC=clang
ubsan: CFLAGS += -fsanitize=address,undefined
ubsan: LDFLAGS += -fsanitize=address,undefined
ubsan: clean benchmark t/00util/api_regress examples/tsnc_sync examples/tsnc_multipath examples/blockchain_gossip
	./benchmark 1280 50000 0
	./t/00util/api_regress
	./examples/tsnc_sync
	./examples/tsnc_multipath
	./examples/blockchain_gossip

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
	$(MAKE) libnanorq_core.a CPPFLAGS="$(CPPFLAGS) -DNANORQ_NO_LIBC"
	@echo "--- Undefined symbols in libnanorq_core.a ---"
	@if nm -u libnanorq_core.a | grep -E '\b(malloc|calloc|realloc|free|posix_memalign|__assert_fail)\b'; then \
		echo "FAIL: libc symbols found in embedded build"; exit 1; \
	else \
		echo "PASS: no libc allocator/assert symbols found"; \
	fi

valgrind: CPPFLAGS=-Wall -Iinclude -Ideps/ -fPIC
valgrind: CFLAGS = -O0 -g -std=c11
valgrind: clean $(TEST_UTILS) $(EXAMPLES)
	valgrind --leak-check=full --errors-for-leak-kinds=definite --error-exitcode=2 ./t/00util/api_regress
	valgrind --leak-check=full --errors-for-leak-kinds=definite --error-exitcode=2 ./examples/tsnc_sync > /dev/null
	valgrind --leak-check=full --errors-for-leak-kinds=definite --error-exitcode=2 ./examples/tsnc_multipath > /dev/null
	valgrind --leak-check=full --errors-for-leak-kinds=definite --error-exitcode=2 ./examples/blockchain_gossip > /dev/null
	valgrind --error-exitcode=2 ./t/00util/hdpcgen  500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/matgen   500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/precond  500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/repgen   500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/ult      500 > /dev/null
	valgrind --error-exitcode=2 ./t/00util/schedgen 500 > /dev/null
	valgrind --error-exitcode=2 ./examples/encode   500 64 10 t/assets/sample_data/raw > /dev/null

check-interop:
	./t/interop_harness.sh
