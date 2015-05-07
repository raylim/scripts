#!/usr/bin/env perl

use strict;
use warnings;

use Cwd qw/ realpath /;
use IO::Socket qw/ AF_UNIX SOCK_STREAM /;
use Path::Class qw/ file /;
use File::Temp();

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
my $scriptFile = File::Temp->new(TEMPLATE => 'tempXXXXX', DIR => '/home/limr/share/tmp', SUFFIX => '.sge');
while (<STDIN>) {
    print $scriptFile $_;
}
close $scriptFile;

my $socketPath = file(realpath($0))->dir()->file('socket');
my $client = IO::Socket->new(
    Domain => AF_UNIX,
    Type => SOCK_STREAM,
    Peer => $socketPath,
) or die("Can't connect to server socket $!\n");

print $client $args;
print $client getcwd();
print $client $scriptFile->filename;

my $exitCode;
while (<$client>) {
    if (/^Error:/) {
        print STDERR;
        $exitCode = -1;
    } elsif (/^Code:/) {
        ($exitCode) = $_ =~ /Exit: (\d+)/;
    }
}

if ($opt{o} && (!-e $opt{o} || !-s $opt{o})) {
    sleep 60; # wait for file system to update
    system("rm $opt{o}");
    #print "File not removed\n" if (-e $opt{o});
    die "$opt{o}: file is size 0";
}

exit $exitCode;
