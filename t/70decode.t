use strict;
use warnings;
use Test::More;
use t::Util;

sub run_decode {
    my ($K, $drops) = @_;
    return run_prog("./examples/decode $K 64 $drops @{[ASSETS_DIR]}/sample_data/raw");
}

subtest "schedules" => sub {
    my @ks = (10, 50, 100, 500, 1000, 5000, 10000, 50000, 56403);
    foreach (@ks) {
        my $drops = int($_ * 0.1);
        $drops = 1 if $drops == 0;
        my ($timing, $resp) = run_decode($_, $drops);
        ok $? == 0, "decode K: $_ drops: $drops";
        if ($? != 0) {
            diag("Error output: $timing");
        }
    }
};

done_testing();
