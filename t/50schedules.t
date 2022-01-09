use strict;
use warnings;
use Test::More;
use t::Util;
use File::Temp qw(tempfile tempdir);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

sub run_schedgen {
    my ($K) = @_;
    return run_prog("./t/00util/schedgen $K");
}

subtest "schedules" => sub {
    my @ks = (10, 50, 100, 500, 1000);
    foreach (@ks) {
        my $expected = md5_file("@{[ASSETS_DIR]}/schedules/K_$_.txt");
        my $resp = run_schedgen($_);
        is md5_hex($resp), $expected, "schedule K: $_";
    }
};
        
done_testing();
