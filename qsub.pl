#!/usr/bin/env perl
# qsub wrapper script

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('hd:', \%opt);

my $usage = <<ENDL;
Usage: perl qsub.pl -d [sge dir] [qsub cmd]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};


my $qsub = shift @ARGV; 

my $args = join " ", @ARGV;

my $qsubCmd = "$qsub $args";
my $qsubOut = qx($qsubCmd);

$qsubOut =~ m/Your job (\d+) /;
my $jobId = $1;
print "$qsubOut";

do {
    my $qstat = qx($opt{d}/qstat -j $jobId 2>&1);
    sleep 30;
} until ($? != 0);

my $exitStatus = qx($opt{d}/qacct -j $jobId | grep exit_status);
#print $exitStatus . "\n";
my $exit = 1;
if ($exitStatus =~ m/exit_status\s+(\d+)/) {
    $exit = $1;
}
#print "exit code: $exit\n";
exit $exit;

