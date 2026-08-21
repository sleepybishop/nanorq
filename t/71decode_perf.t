use strict;
use warnings;
use Test::More;
use t::Util;

sub run_decode {
    my($K, $drops) = @_;
    return run_prog("./examples/decode $K 1280 $drops /dev/urandom");
}

subtest "schedules" => sub {
    my $min_mbps = $ENV{NANORQ_MIN_MBPS} // 1;
    my @ks = (10, 50, 100, 500, 1000, 5000, 10000, 50000, 56403);
    foreach (@ks) {
        my $drops = int($_ * 0.1);
        $drops = 1 if $drops == 0;
        my($timing, $resp) = run_decode($_, $drops);
        ok $? == 0, "decode exited cleanly for K: $_";
        my ($calc, $ops, $bytes) = $timing =~ /calc: ([0-9\.]+)s ops: ([0-9\.]+)s bytes: ([0-9]+)/;
        if (defined $calc && defined $ops && defined $bytes) {
            ok $ops > 0 && ($calc + $ops) > 0 && $bytes > 0,
                "timing values are positive for K: $_" or next;
            my $precalc_mbps = int((8 * $bytes) / ($ops * 1000000));
            my $decode_mbps = int((8 * $bytes) / (($calc + $ops) * 1000000));
            cmp_ok $precalc_mbps, '>=', $min_mbps,
                "PRECALC K: $_ T: 1280 Mbps: $precalc_mbps";
            cmp_ok $decode_mbps, '>=', $min_mbps,
                "DECODE K: $_ T: 1280 Mbps: $decode_mbps";
        } else {
            ok 0, "Failed to parse timing for K: $_";
            diag("Error output: $timing");
        }
    }
};

done_testing();
