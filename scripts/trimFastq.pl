#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $stdout = 0;
my $trim5 = 0;
my $trim3 = 0;
my $FileRead_In;
my $FileRead_Out;

my $usage = "
$0  --trim5 INT --trim3 INT [ - | interleaved-infile  interleaved-outfile ]

each files can be FastQ or gzipped FastQ (*.gz) format

-             read input from STDIN, write output to STDOUT

--trim5 INT   trim INT bases from the 5' end of each read [default $trim5]
--trim3 INT   trim INT bases from the 3' end of each read [default $trim3]

";

GetOptions(""        => \$stdout,
           "trim5=d" => \$trim5,
           "trim3=d" => \$trim3) or die "$usage";

die "must specify input and output or - for stdin/stdout\n\n$usage" if ! $stdout or (!$ARGV[0] and !$ARGV[1]);
die "must specify --trim5 or --trim3\n\n$usage" if ! $trim5 and ! $trim3;
die "must specify positive --trim5 or --trim3\n\n$usage" if $trim5 < 0 or $trim3 < 0;

$FileRead_In = $ARGV[0];
$FileRead_Out = $ARGV[1];
if ($stdout) {
    *INFILE = *STDIN;
} else {
    $FileRead_In = $ARGV[0];
    open(INFILE, "gzip -f -c -d ${FileRead_In} |") or die "couldn't open input $FileRead_In";
}
if ($stdout) {
    *OUTFILE = *STDOUT;
} else {
    $FileRead_Out = $ARGV[1];
    if ($FileRead_Out =~ /\.gz$/) {
    	open(OUTFILE, "| gzip -c > ${FileRead_Out}") or die "couldn't open output $FileRead_Out";
    } else {
    	open(OUTFILE, "> $FileRead_Out") or die "couldn't open output $FileRead_Out";
    }
}

while(<INFILE>) {
    my $r_1 = $_;
    my $r_2 = <INFILE>;
    my $r_3 = <INFILE>;
    my $r_4 = <INFILE>;
    if ($trim5) {
        $r_2 = substr($r_2, $trim5);
        $r_4 = substr($r_4, $trim5);
    }
    if ($trim3) {
        $r_2 = substr($r_2, 0, -$trim3);
        $r_4 = substr($r_4, 0, -$trim3);
    }
    print OUTFILE $r_1;
    print OUTFILE $r_2;
    print OUTFILE $r_3;
    print OUTFILE $r_4;
}
close(INFILE);
close(OUTFILE);

