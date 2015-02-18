#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $minlen = 0;
my $single = "";
my $stdin = 0;
my ($filename1, $filename2, $filenameIn);

my $usage = "
$0  [ interleaved-file ]  read1-outfile  read2-outfile

all files are FASTQ format

-minlen int   Minimum length of sequence to keep.  Both members of a pair 
              are dropped if one is too short [default: $minlen]
-single file  Store reads which meet the -minlen criterion but have a mate
              that did not meet it into this file.  The mate that was too 
              short is dropped.
-             Read input from STDIN.
";

GetOptions("minlen=i" => \$minlen, 
           "single=s" => \$single,
           "" => \$stdin) or die "$usage";

die "$usage" if !$ARGV[0] and !$ARGV[1];

if ($stdin) {
    *INFILE = *STDIN;
    $filename1 = $ARGV[0];
    $filename2 = $ARGV[1];
} else {
    $filenameIn = $ARGV[0];
    $filename1 = $ARGV[1];
    $filename2 = $ARGV[2];
    #if ($filenameIn =~ /\.gz$/) {  # unnecessary to check with these gzip arguments
    	open(INFILE, "gzip -f -c -d ${filenameIn} |") or die "couldn't open $filenameIn";
    #} else {
    #	open(INFILE, "< $filenameIn") or die "couldn't open $filenameIn";
    #}
}
if ($filename1 =~ /\.gz$/) {
    open(FILE1, "| gzip -c > $filename1") or die "couldn't open output $filename1";
} else {
    open(FILE1, "> $filename1") or die "couldn't open output $filename1";
}
if ($filename2 =~ /\.gz$/) {
    open(FILE2, "| gzip -c > $filename2") or die "couldn't open output $filename2";
} else {
    open(FILE2, "> $filename2") or die "couldn't open output $filename2";
}

if ($single) {
    if ($single =~ /\.gz$/) {
        open(SINGLEFILE, "| gzip -c > $single") or die "couldn't open output $single";
    } else {
        open(SINGLEFILE, "> $single") or die "couldn't open output $single";
    }
}

while(<INFILE>) {
        my $f1l1 = $_;
        my $f1l2 = <INFILE>;
        my $f1l3 = <INFILE>;
        my $f1l4 = <INFILE>;
        my $f2l1 = <INFILE>;
        my $f2l2 = <INFILE>;
        my $f2l3 = <INFILE>;
        my $f2l4 = <INFILE>;

        my $s1 = $f1l2;
        chomp $s1;
        my $s2 = $f2l2;
        chomp $s2;
        if (length($s1) >= $minlen && length($s2) >= $minlen) {
		print FILE1 $f1l1;
		print FILE1 $f1l2;
		print FILE1 $f1l3;
		print FILE1 $f1l4;
		print FILE2 $f2l1;
		print FILE2 $f2l2;
		print FILE2 $f2l3;
		print FILE2 $f2l4;
	} elsif ($single) {
		if (length($s1) >= $minlen) {
			print SINGLEFILE $f1l1;
			print SINGLEFILE $f1l2;
			print SINGLEFILE $f1l3;
			print SINGLEFILE $f1l4;
		} elsif (length($s2) >= $minlen) {
			print SINGLEFILE $f2l1;
			print SINGLEFILE $f2l2;
			print SINGLEFILE $f2l3;
			print SINGLEFILE $f2l4;
		}
	}
}

