#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ---------------------------------------------------------------------------
# Section 1 – grid tests: varied (file size, packet size, loss rate)
# ---------------------------------------------------------------------------
SIZES=(
    1           # 1 byte (edge case)
    16          # extremely small
    1300        # single packet slightly > packet size
    65536       # 64 KB, multiple of mmap size
    131072      # 128 KB
    1000000     # ~1 MB
)

PACKET_SIZES=(
    16          # extremely small packet
    500         # normal MTU-ish
    1280        # standard test
    1400        # near full MTU
)

LOSS_RATES=(
    0.0         # perfect network
    1.0         # minimal loss
    10.0        # moderate loss
    50.0        # brutal network conditions
    99.0        # nearly everything dropped (will generate ton of repair packets)
)

# ---------------------------------------------------------------------------
# Section 2 – K-targeted tests: F = K × T to hit exact source-symbol counts
#
# T=1280 is used so K comes out exactly.  Large-K cases (≥5000) only run at
# 0% and 10% loss to keep wall-clock time reasonable.
# ---------------------------------------------------------------------------
# K values to exercise: 10 50 100 500 1000 5000 10000 50000
K_TARGETS=(10 50 100 500 1000 5000 10000 50000)
K_PACKET_SIZE=1280   # T used for K-targeted tests
K_LOSS_RATES_SMALL=(0.0 1.0 10.0 50.0 99.0)   # K ≤ 1000
K_LOSS_RATES_LARGE=(0.0 10.0)                  # K ≥ 5000 (fewer rates for speed)

# ---------------------------------------------------------------------------
echo "Building nanorq..."
(cd .. && make) > /dev/null

echo "Building raptorq_interop..."
cd raptorq_interop
[ -f ~/.cargo/env ] && source ~/.cargo/env || true
cargo build --release --quiet
cd ..

RAPTORQ_INTEROP="raptorq_interop/target/release/raptorq_interop"

if [ ! -f "$RAPTORQ_INTEROP" ]; then
    echo "ERROR: raptorq_interop binary not found at $RAPTORQ_INTEROP"
    exit 1
fi

mkdir -p interop_test_dir
cd interop_test_dir

PASS=0
FAIL=0

# Helper: run one bidirectional test case; updates PASS/FAIL in caller scope
run_test() {
    local size=$1 psize=$2 loss=$3 label=$4

    head -c "$size" /dev/urandom > test_in.bin

    # 1. nanorq encode -> raptorq decode
    if ! ../../encode test_in.bin "$psize" "$loss" 5 > encode.log 2>&1; then
        echo "FAILED [nanorq encode] $label"
        FAIL=$((FAIL + 1)); rm -f data.rq test_out.bin encode.log decode.log test_in.bin; return
    fi
    if ! ../"$RAPTORQ_INTEROP" decode data.rq test_out.bin > decode.log 2>&1; then
        echo "FAILED [raptorq decode] $label"
        FAIL=$((FAIL + 1)); rm -f data.rq test_out.bin encode.log decode.log test_in.bin; return
    fi
    if ! cmp -s test_in.bin test_out.bin; then
        echo "FAILED nanorq->raptorq: $label"
        FAIL=$((FAIL + 1)); rm -f data.rq test_out.bin encode.log decode.log test_in.bin; return
    fi
    rm -f data.rq test_out.bin encode.log decode.log
    PASS=$((PASS + 1))

    # 2. raptorq encode -> nanorq decode
    if ! ../"$RAPTORQ_INTEROP" encode test_in.bin data.rq "$psize" "$loss" 5 > encode.log 2>&1; then
        echo "FAILED [raptorq encode] $label"
        FAIL=$((FAIL + 1)); rm -f data.rq test_out.bin encode.log decode.log test_in.bin; return
    fi
    if ! ../../decode test_out.bin > decode.log 2>&1; then
        echo "FAILED [nanorq decode] $label"
        FAIL=$((FAIL + 1)); rm -f data.rq test_out.bin encode.log decode.log test_in.bin; return
    fi
    if ! cmp -s test_in.bin test_out.bin; then
        echo "FAILED raptorq->nanorq: $label"
        FAIL=$((FAIL + 1)); rm -f data.rq test_out.bin encode.log decode.log test_in.bin; return
    fi
    rm -f data.rq test_out.bin encode.log decode.log test_in.bin
    PASS=$((PASS + 1))
}

# ---------------------------------------------------------------------------
echo ""
echo "=== Section 1: grid tests ==="
for size in "${SIZES[@]}"; do
    for psize in "${PACKET_SIZES[@]}"; do
        for loss in "${LOSS_RATES[@]}"; do
            echo "Testing size: $size bytes, packet size: $psize bytes, loss: $loss%"
            run_test "$size" "$psize" "$loss" "size $size, psize $psize, loss $loss%"
        done
    done
done

# ---------------------------------------------------------------------------
echo ""
echo "=== Section 2: K-targeted tests (T=${K_PACKET_SIZE}) ==="
for k in "${K_TARGETS[@]}"; do
    fsize=$(( k * K_PACKET_SIZE ))
    if (( k >= 5000 )); then
        rates=("${K_LOSS_RATES_LARGE[@]}")
    else
        rates=("${K_LOSS_RATES_SMALL[@]}")
    fi
    for loss in "${rates[@]}"; do
        echo "Testing K=$k (size: $fsize bytes, packet size: ${K_PACKET_SIZE} bytes, loss: $loss%)"
        run_test "$fsize" "$K_PACKET_SIZE" "$loss" "K=$k, loss $loss%"
    done
done

# ---------------------------------------------------------------------------
echo ""
echo "Results: $PASS passed, $FAIL failed"

if [ $FAIL -gt 0 ]; then
    exit 1
fi

echo "All interoperability tests passed!"
