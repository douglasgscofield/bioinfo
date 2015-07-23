#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $stdout = 0;
my $subset = 0;
my $FileRead_1;
my $FileRead_2;
my $FileRead_Out;

my $usage = "
$0  [ - ]  read1-infile  read2-infile [ interleaved-outfile ]

each files can be FastQ or gzipped FastQ (*.gz) format

-             write output to STDOUT
";

GetOptions("subset=i" => \$subset,
           ""         => \$stdout)
or die "$usage";

die "$usage" if !$ARGV[0] and !$ARGV[1];

$FileRead_1 = $ARGV[0];
$FileRead_2 = $ARGV[1];
if ($stdout) {
    *OUTFILE = *STDOUT;
} else {
    $FileRead_Out = $ARGV[2];
    if ($FileRead_Out =~ /\.gz$/) {
    	open(OUTFILE, "| gzip -c > ${FileRead_Out}") or die "couldn't open output $FileRead_Out";
    } else {
    	open(OUTFILE, "> $FileRead_Out") or die "couldn't open output $FileRead_Out";
    }
}
open(FILE1, "gzip -f -c -d ${FileRead_1} |") or die "couldn't open input $FileRead_1";
open(FILE2, "gzip -f -c -d ${FileRead_2} |") or die "couldn't open input $FileRead_2";

while(<FILE1>) {
	print OUTFILE $_;
	$_ = <FILE1>;
	print OUTFILE $_; 
	$_ = <FILE1>;
	print OUTFILE $_; 
	$_ = <FILE1>;
	print OUTFILE $_; 

	$_ = <FILE2>;
	print OUTFILE $_; 
	$_ = <FILE2>;
	print OUTFILE $_;
	$_ = <FILE2>;
	print OUTFILE $_;
	$_ = <FILE2>;
	print OUTFILE $_;

        if ($subset) {
		--$subset;
		print STDERR "subset $subset\n" if ($subset % 100000 == 0);
		last if $subset <= 0;
	}
}
close(FILE1);
close(FILE2);
close(OUTFILE);

