OBJ=\
deps/obl/oblas_lite.o\
lib/params.o\
lib/chooser.o\
lib/precode.o\
lib/rand.o\
lib/tuple.o\
lib/uvec.o\
lib/nanorq.o\
lib/ops.o

TEST_UTILS=\
t/00util/matgen\
t/00util/repgen\
t/00util/hdpcgen\
t/00util/precond\
t/00util/ult\
t/00util/schedgen\
t/00util/bounds

EXAMPLES=\
examples/encode\
examples/decode

CFLAGS   = -O3 -g -std=c11 -Wall -Iinclude -Ideps/ -fPIC -DNDEBUG
CFLAGS  += -march=native -funroll-loops -ftree-vectorize -fno-inline -Wno-unused -Wno-sequence-point -fstack-protector-all

all: libnanorq.a $(EXAMPLES)

t/00util/matgen: t/00util/matgen.o $(OBJ)

t/00util/repgen: t/00util/repgen.o $(OBJ)

t/00util/hdpcgen: t/00util/hdpcgen.o $(OBJ)

t/00util/precond: t/00util/precond.o $(OBJ)

t/00util/ult: t/00util/ult.o $(OBJ)

t/00util/schedgen: t/00util/schedgen.o $(OBJ)

t/00util/bounds: t/00util/bounds.o $(OBJ)

examples/encode: CPPFLAGS += -D_DEFAULT_SOURCE
examples/encode: examples/encode.o $(OBJ)

examples/decode: CPPFLAGS += -D_DEFAULT_SOURCE
examples/decode: examples/decode.o $(OBJ)


check: CPPFLAGS=
check: clean $(TEST_UTILS) $(EXAMPLES)
	prove -I. -v t/*.t

check-nolibc: clean
	$(MAKE) libnanorq.a CPPFLAGS="$(CPPFLAGS) -DNANORQ_NO_LIBC"
	$(MAKE) $(TEST_UTILS)
	prove -I. -v t/10pcmat.t t/20repmat.t t/30precond.t t/35ult.t t/40hdpc.t t/50schedules.t

libnanorq.a:
libnanorq.a: $(OBJ) 
	$(AR) rcs $@ $(OBJ) 

clean:
	$(RM) *.gperf *.prof $(TEST_UTILS) $(EXAMPLES) $(OBJ)
	find -name '*.[a,o]' | xargs $(RM)

indent:
	find -name '*.[h,c]' | xargs clang-format -i

scan:
	scan-build --status-bugs $(MAKE) CPPFLAGS=-D_DEFAULT_SOURCE clean $(OBJ) $(TEST_UTILS) $(EXAMPLES)

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

gperf: LDLIBS = -lprofiler -ltcmalloc
gperf: clean ./examples/encode
	CPUPROFILE_FREQUENCY=100000 CPUPROFILE=gperf.prof ./examples/encode 56403 1280 10 /dev/zero > /dev/null
	pprof ./examples/encode gperf.prof --callgrind > callgrind.gperf
	gprof2dot --format=callgrind callgrind.gperf -z main | dot -T svg > gperf.svg

ubsan: CC=clang
ubsan: CFLAGS += -fsanitize=address,undefined
ubsan: LDFLAGS += -fsanitize=address,undefined
ubsan: LDLIBS += -lubsan
ubsan: clean ./examples/encode
	./examples/encode 56403 1280 10 /dev/zero > /dev/null

benchmark: benchmark.o $(OBJ)

bench: benchmark
	@echo "K       encode   precalc  decode  decode-oh5"
	@./benchmark 1280  100 5.0
	@./benchmark 1280  500 5.0
	@./benchmark 1280 1000 5.0
	@./benchmark 1280 5000 5.0
	@./benchmark 1280 10000 5.0
	@./benchmark 1280 50000 5.0

check-embedded:
	$(MAKE) clean
	$(MAKE) libnanorq.a CPPFLAGS="$(CPPFLAGS) -DNANORQ_NO_LIBC"
	@echo "--- Undefined symbols in libnanorq.a ---"
	@nm -u libnanorq.a | grep -E '\b(malloc|calloc|realloc|free|posix_memalign|__assert_fail)\b' && \
		(echo "FAIL: libc symbols found in embedded build" && exit 1) || \
		echo "PASS: no libc allocator/assert symbols found"


