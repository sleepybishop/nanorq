use strict;
use warnings;
use Test::More;
use t::Util;

my ($stderr, $stdout) = run_prog("./t/00util/api_regress");
is $?, 0, "API regression executable exits successfully";
like $stdout, qr/all API regressions passed/, "API regressions completed";
is $stderr, "", "API regressions have no diagnostics";

done_testing();
