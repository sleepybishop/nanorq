OBJ=\
deps/obl/oblas_lite.o\
lib/params.o\
lib/chooser.o\
lib/precode.o\
lib/rand.o\
lib/tuple.o\
lib/uvec.o\
lib/nanorq.o

TEST_UTILS=\
t/00util/matgen\
t/00util/repgen\
t/00util/hdpcgen\
t/00util/precond\
t/00util/schedgen

EXAMPLES=\
examples/encode\
examples/decode

CPPFLAGS := -DOBLAS_AVX2 -DALIGNSZ=32
CFLAGS   = -O3 -g -std=c11 -Wall -Iinclude -Ideps/ -fPIC
CFLAGS  += -march=native -funroll-loops -ftree-vectorize -Wno-unused -Wno-sequence-point

all: libnanorq.a $(EXAMPLES)

t/00util/matgen: t/00util/matgen.o $(OBJ)

t/00util/repgen: t/00util/repgen.o $(OBJ)

t/00util/hdpcgen: t/00util/hdpcgen.o $(OBJ)

t/00util/precond: t/00util/precond.o $(OBJ)

t/00util/schedgen: t/00util/schedgen.o $(OBJ)

examples/encode: examples/encode.o examples/operations.o $(OBJ)

examples/decode: examples/decode.o examples/operations.o $(OBJ)

check: CPPFLAGS=
check: clean $(TEST_UTILS)
	prove -I. -v t/*.t

libnanorq.a:
libnanorq.a: $(OBJ) 
	$(AR) rcs $@ $(OBJ) 

clean:
	$(RM) *.gperf *.prof $(TEST_UTILS) $(EXAMPLES) $(OBJ)
	find -name '*.[a,o]' | xargs $(RM)

indent:
	find -name '*.[h,c]' | xargs clang-format -i

scan:
	scan-build $(MAKE) CPPFLAGS= clean $(OBJ) $(TEST_UTILS) $(EXAMPLES)

gperf: LDLIBS = -lprofiler -ltcmalloc
gperf: clean ./examples/encode
	CPUPROFILE_FREQUENCY=100000 CPUPROFILE=gperf.prof ./examples/encode 56403 1280 10 /dev/zero > /dev/null
	pprof ./examples/encode gperf.prof --callgrind > callgrind.gperf
	gprof2dot --format=callgrind callgrind.gperf -z main | dot -T svg > gperf.svg

ubsan: CC=clang
ubsan: CFLAGS += -fsanitize=undefined
ubsan: LDLIBS += -lubsan
ubsan: clean ./examples/encode
	./examples/encode 56403 1280 10 /dev/zero > /dev/null

