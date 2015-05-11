#!/usr/bin/env perl

use strict;
use warnings;

use Cwd qw/ realpath /;
#use IO::Socket qw/ AF_UNIX SOCK_STREAM /;
use IO::Socket::INET;
use Path::Class qw/ file /;
use File::Temp();
use Cwd;
use Storable qw/ nfreeze /;

use Getopt::Std;
my %opt;
getopts('ho:', \%opt);

my $usage = <<ENDL;
Usage: perl qsubClient.pl -h -- [qsub args]
    -o [file]: check file for non-zero size
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my $args = join " ", @ARGV;
my $socketPath = file($opt{s});
#my $client = IO::Socket->new(
#Domain => AF_UNIX,
#Type => SOCK_STREAM,
#Peer => $socketPath,
#Timeout => 30,
#) or die("Can't connect to server socket: $!\n");
#

my $client = IO::Socket::INET->new(
    PeerHost => 'localhost',
    PeerPort => '34383',
    Proto => 'tcp',
) or die "Can't connect to server socket: $!\n";
eval 'END { close $client } 1' or die $@;


#print "Connected to server $socketPath\n";
#print "Sending server args: $args\n";
print $client $args . "\n";
#print "Sending server cwd: " . getcwd() . "\n";
print $client getcwd() . "\n";

my $scriptFile = File::Temp->new(TEMPLATE => 'tempXXXXX', DIR => '/home/limr/share/tmp', SUFFIX => '.sge', UNLINK => 0);
while (my $line = <STDIN>) {
    print $scriptFile $line;
}
close $scriptFile;
#print "Sending server script: " . $scriptFile->filename . "\n";
print $client $scriptFile->filename . "\n";

my $exitCode = -1;
while (<$client>) {
    if (/^Error:/) {
        print STDERR;
        $exitCode = -1;
        last;
    } elsif (/^Code:/) {
        ($exitCode) = $_ =~ /Code: (\d+)/;
        last;
    }
}

if ($opt{o} && (!-e $opt{o} || !-s $opt{o})) {
    sleep 60; # wait for file system to update
    system("rm $opt{o}");
    #print "File not removed\n" if (-e $opt{o});
    die "$opt{o}: file is size 0";
}

exit $exitCode;
