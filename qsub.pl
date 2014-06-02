#!/usr/bin/env perl
# qsub wrapper script

use strict;
use warnings;

use Schedule::DRMAAc qw/ :all /;
use File::Temp();
use Cwd;

use Getopt::Std;
my %opt;
getopts('h', \%opt);

my $usage = <<ENDL;
Usage: perl qsub.pl -h -- [qsub args]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my $scriptFile = File::Temp->new(TEMPLATE => 'tempXXXXX', DIR => '/tmp', SUFFIX => '.sge');

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
    die "Received interrupt: terminating job\n";
}

$SIG{INT} = \&signalHandler;
$SIG{TERM} = \&signalHandler;

# loop to give a chance to receive sigint/sigterms
my $stat;
do {
    ($error, my $jobidOut, $stat, $diagnosis) = drmaa_wait($jobid, 10);
} until ($error != $DRMAA_ERRNO_EXIT_TIMEOUT);

# pull all exit-related codes
($error, my $exitStatus, $diagnosis) = drmaa_wexitstatus($stat);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
($error, my $aborted, $diagnosis) = drmaa_wifaborted( $stat );
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
($error, my $signaled, $diagnosis ) = drmaa_wifsignaled( $stat );
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
($error, my $coreDumped, $diagnosis ) = drmaa_wcoredump( $stat );
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_exit();
die drmaa_strerror($error) . "\n" . $diagnosis if $error;


exit $exitStatus + $aborted + $signaled + $coreDumped;
