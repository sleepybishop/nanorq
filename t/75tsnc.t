use strict;
use warnings;
use Test::More;
use t::Util;

for my $demo (qw(tsnc_sync tsnc_multipath blockchain_gossip)) {
    my ($stderr, $stdout) = run_prog("./examples/$demo");
    ok $? == 0, "$demo exits successfully";
    like $stdout, qr/\[Verification\] SUCCESS:/,
        "$demo verifies decoded output";
}

done_testing();
