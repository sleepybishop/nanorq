# nanorq
[![CI](https://github.com/sleepybishop/nanorq/actions/workflows/ci.yml/badge.svg?branch=stable)](https://github.com/sleepybishop/nanorq/actions/workflows/ci.yml)

nanorq is a compact, performant implementation of the raptorq fountain code capable of reaching multi-gigabit speeds on a single core.

## API
  - High-level API `nanorq.h` abstracts allocation, alignment, and I/O handling from the user.
  - Low-level API `nanorq_core.h` provides more flexibility and can be used without any libc dependency[^1], but requires more integration from the user.

## Performance[^2]
![](graph.png)

## Use cases
  - firmware deployment / software updates
  - video streaming
  - large data transfers across high latency links

[^1]: When compiled with `NDEBUG`. Test suite requires stdlib.
[^2]: Build autodetects and optimizes for the host CPU architecture (using -march=native).
