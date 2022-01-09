use strict;
use warnings;
use Test::More;
use t::Util;
use File::Temp qw(tempfile tempdir);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

sub run_hdpcgen {
    my ($K) = @_;
    return run_prog("./t/00util/hdpcgen $K");
}

subtest "hdpc" => sub {
    my @ks = (10, 50, 100, 500, 1000, 5000, 10000, 56403);
    foreach (@ks) {
        my $resp = run_hdpcgen($_);
        ok $resp =~ /^HDPC\[(\d+)x(\d+)\]/, "matches header $1x$2";
        my $expected = md5_file("@{[ASSETS_DIR]}/mats/hdpc/K_$_.txt");
        is md5_hex($resp), $expected, "hdpc K: $_";
    }
};
        
done_testing();
