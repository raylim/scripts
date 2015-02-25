#!/usr/bin/env perl
# usage: perl checkMD5.pl {file}

use strict;

my $passes;
my $fails;
foreach (@ARGV) {
    my $md5file = $_ . ".md5";
    if (! -f $md5file) {
        die "$md5file does not exist\n";
    }
    chomp(my $md5 = `cut -f1 -d' ' $md5file`);
    chomp(my $newMd5 = `md5 $_ | sed 's/.*= //'`);
    if ($md5 eq $newMd5) {
        print "PASS: $_ ($md5)\n";
        $passes++;
    } else {
        print "FAIL: $_\n\t$md5 != $newMd5\n";
        $fails++;
    }
}
if ($fails) {
    print "$passes / " . scalar @ARGV . " PASSED\n";
    exit 1;
} else {
    print "ALL PASSED\n";
}

