use strict;
use warnings;
use Test::More;
use t::Util;
use File::Temp qw(tempfile tempdir);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

sub run_precond {
    my ($K) = @_;
    return run_prog("./t/00util/precond $K");
}

subtest "precond" => sub {
    my @ks = (10, 50, 100, 500, 1000);
    foreach (@ks) {
        my $resp = run_precond($_);
        ok $resp =~ /^A\[(\d+)x(\d+)\]/, "matches header $1x$2";
        my $expected = md5_file("@{[ASSETS_DIR]}/mats/apc/K_$_.txt");
        is md5_hex($resp), $expected, "precond K: $_";
    }
};
        
done_testing();
