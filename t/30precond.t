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

sub check_lt {
    my @lines = reverse(split("\n", shift @_));
    my ($nz, $last, $row) = (0,0,0);
    my $single = shift @lines;
    if ($single !~ /^(0+)1$/) {
        return 0;
    }
    $nz = length($1);
    foreach (@lines) {
        if ($_ =~ /^(0+)/) {
            $row++;
            $nz--;
            last if $nz == 0;
            if ($nz != length($1)) {
                diag Dumper($nz, length($1));
                return 0;
            }
        }
    }
    return $row;
}

subtest "precond" => sub {
    my @ks = (10, 50, 100, 500, 1000);
    foreach (@ks) {
        my $resp = run_precond($_);
        ok $resp =~ /^U\[(\d+)x(\d+)\]/, "matches header $1x$2";
        my $expected = md5_file("@{[ASSETS_DIR]}/mats/u/K_$_.txt");
        is md5_hex($resp), $expected, "precond K: $_";
        ok check_lt($resp), "lower triangular";
    }
};
        
done_testing();
