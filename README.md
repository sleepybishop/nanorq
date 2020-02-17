# nanorq
this is a c99 port of raptorq codes aka rfc6330 with a focus on performace and code size.

it was originally inspired by libraptorq and velopyraptorq but diverged significantly in it's approach and api.

afaik its the fastest public implementation however performance could still be much improved by optimizing the preconditioner and supporting inverse updates so the schedule could be cached when decoding multiple source blocks.
