use strict;
use warnings;
use Test::More;
use t::Util;
use File::Temp qw(tempfile tempdir);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

sub run_encode
{
    my($K) = @_;
    return run_prog("./examples/encode $K 1280 10 /dev/urandom");
}

subtest "schedules" => sub
{
    my $min_mbps = $ENV{NANORQ_MIN_MBPS} // 1;
    my @ks = (10, 50, 100, 500, 1000, 5000, 10000, 50000, 56403);
    foreach (@ks) {
        my $expected = md5_file("@{[ASSETS_DIR]}/examples/encode/K_$_.txt");
        my($timing, $resp) = run_encode($_);
        ok $? == 0, "encode exited cleanly for K: $_";
        my ($calc, $ops, $bytes) = $timing =~ /calc: ([0-9\.]+)s ops: ([0-9\.]+)s bytes: ([0-9]+)/;
        ok defined($calc) && defined($ops) && defined($bytes),
            "parsed timing for K: $_" or next;
        ok $ops > 0 && ($calc + $ops) > 0 && $bytes > 0,
            "timing values are positive for K: $_" or next;
        my $precalc_mbps = int((8 * $bytes) / ($ops * 1000000));
        my $encode_mbps = int((8 * $bytes) / (($calc + $ops) * 1000000));
        cmp_ok $precalc_mbps, '>=', $min_mbps,
            "PRECALC K: $_ T: 1280 Mbps: $precalc_mbps";
        cmp_ok $encode_mbps, '>=', $min_mbps,
            "ENCODE K: $_ T: 1280 Mbps: $encode_mbps";
    }
};

done_testing();
