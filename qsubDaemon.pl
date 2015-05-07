#!/usr/bin/env perl
# qsub daemon 

use strict;
use warnings;

use Proc::Daemon;
use threads;
use threads::shared;
use Thread qw/ async /;
use Thread::Queue;
use Schedule::DRMAAc qw/ :all /;
use File::Temp();
use IO::Socket qw/ AF_UNIX SOCK_STREAM SOMAXCONN/;
use Path::Class qw/ file /;
use Cwd;

Proc::Daemon::Init;

sub done { exit(0); }
$SIG{TERM} = \&done;
$SIG{INT} = \&done;

# process job
my $jobQueue = new Thread::Queue;
my $thr = async {
    while (my $clientJob = $jobQueue->dequeue) {
        my ($client, $jobid) = ($clientJob->[0], $clientJob->[1]);
        
        # delete jobs of disconnected clients
        unless ($client->connected) {
            my ($error, $diagnosis) = drmaa_control($jobid, $DRMAA_CONTROL_TERMINATE);
            print "Error: " .drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;
        }

        my ($error, $jobidOut, $stat, $diagnosis) = drmaa_wait($jobid, 3);
        if ($error == $DRMAA_ERRNO_EXIT_TIMEOUT) {
            # tell client job is complete
            # pull all exit-related codes
            ($error, my $exitStatus, $diagnosis) = drmaa_wexitstatus($stat);
            print $client "Error: " .drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;
            ($error, my $aborted, $diagnosis) = drmaa_wifaborted( $stat );
            print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;
            ($error, $diagnosis) = drmaa_exit();
            print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;
            }
            print $client "Code: " . $exitStatus + $aborted + $signaled + $coreDumped;
        } else {
            # requeue the job
            $jobQueue->enqueue($clientJob);
        }
    }
};

my $socketPath = file(realpath($0))->dir()->file('socket');
 
{
    my $server = new IO::Socket->new (
        Domain => AF_UNIX,
        Type => SOCK_STREAM,
        Local => $socketPath,
        Listen => SOMAXCONN,
    ) or die "Unable to create server soket: $!\n";
    eval 'END { unlink $socketPath} 1' or die $@;

    my ($error, $diagnosis) = drmaa_init(undef);
    die drmaa_strerror($error) . "\n" . $diagnosis if $error;

    while (my $client = $socket->accept()) {
        my $clientArgs = <$clientsocket>;
        my $clientCwd = <$clientsocket>;
        my $clientScript = <$clientsocket>;

        my $scriptFile = File::Temp->new(TEMPLATE => 'tempXXXXX', DIR => '/home/limr/share/tmp', SUFFIX => '.sge');
        print $scriptFile $clientScript;
        close $scriptFile;
        
        ($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $scriptFile->filename);
        print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and continue if $error;

        ($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, $args);
        print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and continue if $error;

        ($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_WD, getcwd());
        print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and continue if $error;

        ($error, my $jobid, $diagnosis) = drmaa_run_job($jt);
        print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and continue if $error;

        ($error, $diagnosis) = drmaa_delete_job_template($jt)
        print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and continue if $error;

        $jobQueue->enqueue([ $client, $jobid ]);
    }
}
