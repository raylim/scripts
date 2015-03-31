#!/usr/bin/env perl
# usage: perl checkMD5.pl {file}

use strict;

use Getopt::Std;

our($opt_f);
getopts('f:');

my @files;
if (defined $opt_f) {
    open IN, "<$opt_f";
    @files = <IN>;
    close IN;
} else {
    @files = @ARGV;
}

my $passes;
my $fails;
foreach (@files) {
    chomp;
    my $md5file = $_ . ".md5";
    if (! -f $md5file) {
        print "$md5file does not exist\n";
        next;
    }
    print "$_: ";
    chomp(my $md5 = `cut -f1 -d' ' $md5file`);
    chomp(my $newMd5 = `md5 $_ | sed 's/.*= //'`);
    if ($md5 eq $newMd5) {
        print "PASS ($md5)\n";
        $passes++;
    } else {
        print "FAIL\n\t$md5 != $newMd5\n";
        $fails++;
    }
}
if ($fails) {
    print "$passes / " . scalar @ARGV . " PASSED\n";
    exit 1;
} else {
    print "ALL PASSED\n";
}

