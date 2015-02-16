#!/usr/bin/env perl

use Getopt::Long;

my $interleaved = 0;
my $stdin = 0;

my $usage = "
$0  

all files are FASTQ format or gzipped FASTQ (*.gz)

-interleaved  file is interleaved
-             read input from STDIN
";

GetOptions("interleaved" => \$interleaved, 
           "" => \$stdin) or { die "$usage" };

die "$usage" if !$stdin and !$ARGV[0];
my $File = $ARGV[0];

if ($stdin) {
    *INFILE = *STDIN;
} else {
    open(INFILE, "gzip -f -c -d ${File} |") or die "couldn't open $File";
}

my $nseq = 0;
my $nbase = 0;

while(<INFILE>) {
    my $l = $_;  # @name
    my $l2 = <INFILE>;  # sequence
    $l = <INFILE>;  # +
    $l = <INFILE>;  # quality
    ++$nseq;
    $nbase += (length($l2) - 1); # chomp;
}

print "nseq\t", $nseq, "\n";
print "npair\t", ($nseq/2), "\n" if $interleaved;
print "nbase\t", $nbase, "\n";
print "nMbp\t", ($nbase / 1000000), "\n";
print "meanlen\t", ($nbase / $nseq), "\n";
