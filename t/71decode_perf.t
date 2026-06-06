use strict;
use warnings;
use Test::More;
use t::Util;

sub run_decode {
    my($K, $drops) = @_;
    return run_prog("./examples/decode $K 1280 $drops /dev/urandom");
}

subtest "schedules" => sub {
    my @ks = (10, 50, 100, 500, 1000, 5000, 10000, 50000, 56403);
    foreach (@ks) {
        my $drops = int($_ * 0.1);
        $drops = 1 if $drops == 0;
        my($timing, $resp) = run_decode($_, $drops);
        ok $? == 0, "decode exited cleanly for K: $_";
        my ($calc, $ops, $bytes) = $timing =~ /calc: ([0-9\.]+)s ops: ([0-9\.]+)s bytes: ([0-9]+)/;
        if (defined $calc && defined $ops && defined $bytes) {
            my $precalc_mbps = $ops > 0 ? int((8 * $bytes) / ($ops * 1000000)) : "Inf";
            my $decode_mbps = ($calc + $ops) > 0 ? int((8 * $bytes) / (($calc + $ops) * 1000000)) : "Inf";
            ok 1, "PRECALC K: $_ T: 1280 Mbps: $precalc_mbps";
            ok 1, "DECODE K: $_ T: 1280 Mbps: $decode_mbps";
        } else {
            ok 0, "Failed to parse timing for K: $_";
            diag("Error output: $timing");
        }
    }
};

done_testing();
