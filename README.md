# nanorq
this is a c99 port of raptorq codes aka rfc6330, based off libraptorq and others with a focus on performace and code size.

afaik its the fastest public implementation however a couple of large issues remain.

1. the decoder can handle dependent rows a lot better in the case of non zero overhead.

2. performance could still be much improved by utilizing sparse matrices during precoding.
