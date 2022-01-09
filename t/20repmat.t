use strict;
use warnings;
use Test::More;
use t::Util;
use File::Temp qw(tempfile tempdir);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

sub run_repgen {
    my ($K) = @_;
    return run_prog("./t/00util/repgen $K");
}

sub replace_lines {
    my @lines = split("\n", shift @_);
    my ($from, $to) = @_;
    @lines = @lines[0..$from,($to+1)..$#lines];
    return join("\n", @lines);
}

subtest "repmat" => sub {
    my @ks = (10, 50, 100, 500, 1000);
    foreach (@ks) {
        my $resp = run_repgen($_);
        ok $resp =~ /^A\[(\d+)x(\d+)\],S_H\[(\d+)\|(\d+)\]/, "matches header $1x$2, LDPC: $3, HDPC: $4";
        my $expected = slurp_file("@{[ASSETS_DIR]}/mats/a/K_$_.txt");
        $resp = replace_lines($resp, $3+$4, $3+$4+3);
        $expected = replace_lines($expected, $3+$4, $3+$4+3);
        is $resp, $expected, "repmat K: $_";
    }
};
        
done_testing();
