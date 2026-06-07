use strict;
use warnings;
use Test::More;
use t::Util;

subtest "utilities" => sub {
    my ($stderr, $stdout) = run_prog("./t/00util/test_utils");
    like $stdout, qr/All utility tests passed successfully!/, "utility tests output success message";
    is $stderr, '', "no errors printed to stderr";
};

done_testing();
