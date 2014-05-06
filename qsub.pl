#!/usr/bin/env perl
# qsub wrapper script

use strict;
use warnings;

use Schedule::DRMAAc qw/ :all /;
use File::Temp();
use Cwd;

use Getopt::Std;
my %opt;
getopts('hS:', \%opt);

my $usage = <<ENDL;
Usage: perl qsub.pl -S [shell] -- [qsub args]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my $shell = "/bin/bash";
$shell = $opt{s} if $opt{s};

my $scriptFile = File::Temp->new(TEMPLATE => 'tempXXXXX', SUFFIX => '.sge');

#print $scriptFile->filename . "\n";
print $scriptFile "#!$shell\n";

my $args = join " ", @ARGV;
while (<STDIN>) {
    print $scriptFile $_;
}
close $scriptFile;

my ($error, $diagnosis) = drmaa_init(undef);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, my $jt, $diagnosis) = drmaa_allocate_job_template();
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $scriptFile->filename);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, $args);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_WD, getcwd());
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, my $jobid, $diagnosis) = drmaa_run_job($jt);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

sub signalHandler {
    my ($error, $diagnosis) = drmaa_control($jobid, $DRMAA_CONTROL_TERMINATE);
    die drmaa_strerror($error) . "\n" . $diagnosis if $error;
    die;
}

$SIG{INT} = \&signalHandler;
$SIG{TERM} = \&signalHandler;

($error, my $jobIdOut, my $stat, $diagnosis) = drmaa_wait($jobid, 30);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, my $exited, $diagnosis) = drmaa_wifexited($stat);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

exit $exited;

