#!/usr/bin/env perl

use Getopt::Long;

my $stdin = 0;
my $single = 0;
my $fraction = 0;
my $count = 0;
my $o_count = 0;
my $filenameIn = "";
my $filenameOut = "";
my $n_in = 0;
my $n_out = 0;

my $usage = "
$0  [ fastq-file fastq-file ]

All files are interleaved FASTQ, and may be gzipped (*.gz).  Single-end
reads can be handles with --single

-c|--count INT        number of reads to keep... OR
-f|--fraction FLOAT   fraction of reads to keep; a read pair is
                      selected if a random uniform draw <= FLOAT
-s|--single           reads are single-end
-                     read input from STDIN, write to STDOUT

To handle paired-end reads in separate files for reads 1 and 2, use this
script in a pipe like

    shuffleFastq.pl - read1.fq.gz read2.fq.gz \
    | subsampleReads.pl - -f 0.01 \
    | deshuffleFastq.pl - sub1.fq.gz sub2.fq.gz

";

GetOptions("f|fraction=f" => \$fraction, 
           "c|count=i"    => \$o_count, 
           "s|single"     => \$single, 
           ""             => \$stdin) or { die "$usage" };

die "$usage" if $stdin and ($ARGV[0] or $ARGV[1]);
die "$usage" if ! $stdin and (!$ARGV[0] or !$ARGV[1]);
die "only one of --fraction or --count" if $fraction and $o_count;
die "fraction is strange ($fraction)" if not $o_count and ($fraction <= 0.0 or $fraction >= 1.0);
$count = $o_count;

if ($stdin) {
    *INFILE = *STDIN;
    *OUTFILE = *STDOUT;
} else {
    $filenameIn = $ARGV[0];
    $filenameOut = $ARGV[1];
    if ($filenameIn =~ /\.gz$/) {
        open(INFILE, "gzip -f -c -d ${filenameIn} |") or die "couldn't open $filenameIn";
    } else {
        open(INFILE, "< $filenameIn") or die "couldn't open $filenameIn";
    }
    if ($filenameOut =~ /\.gz$/) {
        open(OUTFILE, "| gzip -c > $filenameOut") or die "couldn't open output $filenameOut";
    } else {
        open(OUTFILE, "> $filenameOut") or die "couldn't open output $filenameOut";
    }
}

my ($f1l1, $f1l2, $f1l3, $f1l4, $f2l1, $f2l2, $f2l3, $f2l4);

while(<INFILE>) {
    $f1l1 = $_;
    $f1l2 = <INFILE>;
    $f1l3 = <INFILE>;
    $f1l4 = <INFILE>;
    if (! $single) {
        $f2l1 = <INFILE>;
        $f2l2 = <INFILE>;
        $f2l3 = <INFILE>;
        $f2l4 = <INFILE>;
    }
    ++$n_in;

    if ($fraction) {
        if (rand() <= $fraction) {
            print OUTFILE $f1l1;
            print OUTFILE $f1l2;
            print OUTFILE $f1l3;
            print OUTFILE $f1l4;
            if (! $single) {
            print OUTFILE $f2l1;
            print OUTFILE $f2l2;
            print OUTFILE $f2l3;
            print OUTFILE $f2l4;
            }
            ++$n_out;
        }
    } elsif ($count) {
        print OUTFILE $f1l1;
        print OUTFILE $f1l2;
        print OUTFILE $f1l3;
        print OUTFILE $f1l4;
        if (! $single) {
        print OUTFILE $f2l1;
        print OUTFILE $f2l2;
        print OUTFILE $f2l3;
        print OUTFILE $f2l4;
        }
        ++$n_out;
        --$count;
        last if not $count;
    }
}

printf STDERR "Input reads: %d  Output reads: %d  Count requested: %d  Fraction requested: %0.5f  Fraction realized: %0.5f\n",
              $n_in, $n_out, $o_count, $fraction, $n_out/$n_in;
