CFLAGS = -O3 -g -std=c11 -Wall -I../include -I../deps -I. -march=native

all: libtsnc.a examples/tsnc_sync examples/tsnc_multipath examples/blockchain_gossip

libtsnc.a: tsnc.o
	ar rcs $@ $^

tsnc.o: tsnc.c tsnc.h
	$(CC) $(CFLAGS) -c $< -o $@

examples/tsnc_sync: examples/tsnc_sync.c libtsnc.a
	mkdir -p examples
	$(CC) $(CFLAGS) $< -o $@ -L. -ltsnc -L../ -lnanorq -lm

examples/tsnc_multipath: examples/tsnc_multipath.c libtsnc.a
	mkdir -p examples
	$(CC) $(CFLAGS) $< -o $@ -L. -ltsnc -L../ -lnanorq -lm

examples/blockchain_gossip: examples/blockchain_gossip.c libtsnc.a
	mkdir -p examples
	$(CC) $(CFLAGS) $< -o $@ -L. -ltsnc -L../ -lnanorq -lm

clean:
	rm -f *.o *.a examples/tsnc_sync examples/tsnc_multipath examples/blockchain_gossip
