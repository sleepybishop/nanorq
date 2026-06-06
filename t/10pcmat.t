use strict;
use warnings;
use Test::More;
use t::Util;
use File::Temp qw(tempfile tempdir);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

sub run_matgen {
    my ($K) = @_;
    return run_prog("./t/00util/matgen $K");
}

subtest "pcmat" => sub {
    my @ks = (10, 50, 100, 500, 1000);
    foreach (@ks) {
        my $resp = run_matgen($_);
        ok $resp =~ /^A\[(\d+)x(\d+)\],S_H\[(\d+)\|(\d+)\]/, "matches header $1x$2, LDPC: $3, HDPC: $4";
        my $expected = md5_file("@{[ASSETS_DIR]}/mats/a/K_$_.txt");
        is md5_hex($resp), $expected, "pcmat K: $_";
    }
};
        
done_testing();
