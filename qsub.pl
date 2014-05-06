#!/usr/bin/env perl
# qsub wrapper script

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('hd:s:', \%opt);

my $usage = <<ENDL;
Usage: perl qsub.pl -s [min sleep seconds] -d [sge dir] [qsub cmd]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};


my $qdelCmd = "qdel";
my $qstatCmd = "qstat";
my $qacctCmd = "qacct";
if ($opt{d}) {
    $qdelCmd = $opt{d} . "/" . $qdelCmd if $opt{d};
    $qstatCmd = $opt{d} . "/" . $qstatCmd if $opt{d};
    $qacctCmd = $opt{d} . "/" . $qacctCmd if $opt{d};
}


my $sleepTime = 5 unless $opt{s};

my $qsub = shift @ARGV; 

my $args = join " ", @ARGV;

my $qsubCmd = "$qsub $args";
my $qsubOut = qx($qsubCmd);

$qsubOut =~ m/Your job (\d+) /;
my $jobId = $1;
print "$qsubOut";


sub signalHandler {
    system "$qdelCmd $jobId";
    die;
}

$SIG{INT} = \&signalHandler;
$SIG{TERM} = \&signalHandler;

do {
    my $qstat = qx($qstatCmd -j $jobId 2>&1);
    sleep $sleepTime + int(rand(10));
} until ($? != 0);

my $exitStatus = qx($qacctCmd -j $jobId | grep exit_status);
#print $exitStatus . "\n";
my $exit = 1;
if ($exitStatus =~ m/exit_status\s+(\d+)/) {
    $exit = $1;
}
#print "exit code: $exit\n";
exit $exit;

