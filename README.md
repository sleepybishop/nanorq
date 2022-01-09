# nanorq
[![CI](https://github.com/sleepybishop/nanorq/actions/workflows/ci.yml/badge.svg?branch=dev)](https://github.com/sleepybishop/nanorq/actions/workflows/ci.yml)

nanorq is a compact, performant implementation of the raptorq fountain code capable of reaching multi-gigabit speeds on a single core.

## Features
  - **No internal memory allocations**. nanorq will ask you to provide the memory it needs.
  - **No stdlib required**. The library itself can be built with `--nostdlib`[^1]
  - **Data privacy**. The library generates operations to run on the data to be coded but the execution of those operations is left to the application.

## Performance[^2]
![](graph.png)

## Use cases
- firmware deployment / software updates
- video streaming
- large data transfers across high latency links

[^1]: When compiled with `NDEBUG`. Test suite requires stdlib.
[^2]: Default build is configured for `AVX2`, adjust Makefile as needed for other archs.
