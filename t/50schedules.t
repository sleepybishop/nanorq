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
    my @lens = (394, 1377, 2259, 11222, 22172);
    foreach (@ks) {
        my $expected = md5_file("@{[ASSETS_DIR]}/schedules/K_$_.txt");
        my $resp = run_schedgen($_);
        my $schedule_len = $resp =~ tr/\n//;
        my $expected_schedule_len = shift @lens;
        is md5_hex($resp), $expected, "schedule K: $_";
        diag $expected_schedule_len;
        is ($schedule_len, $expected_schedule_len, "matched expected length");
    }
};
        
done_testing();
