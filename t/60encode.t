use strict;
use warnings;
use Test::More;
use t::Util;
use File::Temp qw(tempfile tempdir);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

sub run_encode {
    my ($K) = @_;
    return run_prog("./examples/encode $K 64 10 @{[ASSETS_DIR]}/sample_data/raw");
}

subtest "schedules" => sub {
    my @ks = (10, 50, 100, 500, 1000, 5000, 10000, 50000, 56403);
    foreach (@ks) {
        my $expected = md5_file("@{[ASSETS_DIR]}/examples/encode/K_$_.txt");
        my ($timing, $resp) = run_encode($_);
        is md5_hex($resp), $expected, "encode K: $_";
    }
};
        
done_testing();
