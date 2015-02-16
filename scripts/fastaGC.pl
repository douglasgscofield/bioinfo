#!/usr/bin/env perl

###########################################################################
#
# Copyright (c) 2011,2015 Douglas G. Scofield, Uppsala University
#
# Use as you see fit, include the above attribution.  No warranty is implied or
# assumed with this code.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $version = "fastaGC.pl v. 0.1";

=head1 NAME

fastaGC.pl - Compute GC content for various views of FASTA-format sequences

=head1 SYNOPSIS

fastaGC.pl [options] [file ...]

By default fastaGC.pl will print length, GC and base content information across
all FASTA sequences in the input file(s).  More detailed information may be
requested with options.  For all information, GC content is printed in the form
of a tab-separated table suitable for import into whatever other program can
read tables.

Some options can print GC content on a block-by-block basis.  These options
break input sequences into blocks, where each block is one line of FASTA
sequence.  The scale of a block for these options is dependent upon the input
format.  If sequences are broken into 60-character lines, then block size will
be 60 bp.

The options are not mutually exclusive.  All four types of GC content may be
displayed in a single run.

  Options:
    -nototal        do not print information for total input
    -concat         print block-by-block GC content for input sequences
                    as if they were one concatenated sequence
    -block          print mean GC content across all input sequences
                    on a block-by-block basis
    -seq            print GC content for each individual sequence
                    in the input file(s)
    -verbose        print progress (every 1000 sequences), and print
                    headers prior to the output of each table of GC content
    -help, -?       brief help message

=cut

my $printTotalGCStats = 1;
my $printConcatenatedGCStats = 0;
my $printBlockGCStats = 0;
my $printSeqGCStats = 0;
my $verbose = 0;
my $help;

GetOptions(
    "nototal" => sub { $printTotalGCStats = 0; },
    "concat" => \$printConcatenatedGCStats,
    "block" => \$printBlockGCStats,
    "seq" => \$printSeqGCStats,
    "verbose" => \$verbose,
    "help|?" => \$help) or pod2usage(2);

print "GC content of FASTA sequences calculated by $version\n\n" if $help or $verbose;

pod2usage(1) if $help;

my $Nseqs;
my $Nfastalines;
my $Nletters;
my $Nbases;
my $letter;
my $i;
my @ACGTN = (0, 0, 0, 0, 0);
my $seqline = 0;
my $maxseqline = 0;
my $lineGC = 0;
my $lineBases = 0;
my $lineLetters = 0;
my $seqGC = 0;
my $seqBases = 0;
my $seqLetters = 0;
my @seqGC = ();
my @seqBases = ();
my @seqLetters = ();
my $thisseq = "";
my @seqNamesList = ();
my @seqLettersList = ();
my @seqBasesList = ();
my @seqGCContentList = ();

if ($printConcatenatedGCStats) {

	# GC content along the concatenated set of seqs

	print "\n============================================================\n" if $verbose;
	print "GC content along the concatenated set of seqs\n\n" if $verbose;

	print "block\tseq\tGC\tletters\tbases\tGCinclN\tGCexclN\n";
}

while (<>) {
	chomp;
	if (/^>/) {
		if ($thisseq ne "") {
			# update seq stats (seq name, GC content, seq length)
			push @seqNamesList, substr($thisseq, 1);
			push @seqLettersList, $seqLetters;
			push @seqBasesList, $seqBases;
			push @seqGCContentList, $seqGC;
			$seqGC = $seqBases = $seqLetters = 0;
		}
		$thisseq = $_;
		++$Nseqs;
		print "processing seq $Nseqs ...\n" if ($Nseqs % 1000 == 0) and $verbose;
		$maxseqline = $seqline if ($seqline > $maxseqline);
		$seqline = -1;
		next;
	} else {
		++$seqline;
	}
	# print "$seqline\n";
	++$Nfastalines;
	$seqLetters[$seqline] += length;
	$seqBases[$seqline] += length;
	$lineLetters = $lineBases = length;
	$lineGC = 0;
	my @line = split //;
	foreach $letter (@line) {
		if ($letter eq "A") {
			++$ACGTN[0];
		} elsif ($letter eq "C") {
			++$ACGTN[1];
			++$seqGC[$seqline];
			++$lineGC;
		} elsif ($letter eq "G") {
			++$ACGTN[2];
			++$seqGC[$seqline];
			++$lineGC;
		} elsif ($letter eq "T") {
			++$ACGTN[3];
		} else {
			++$ACGTN[4];
			--$seqBases[$seqline];
			--$lineBases;
		}
	}
	$seqGC += $lineGC;
	$seqLetters += $lineLetters;
	$seqBases += $lineBases;

	# these stats are linear across seqs

	if ($printConcatenatedGCStats) {
		if ($lineBases >= 50) {
			printf "%d\t%d\t%d\t%d\t%d\t", 
				$Nfastalines, 
				$Nseqs, 
				$lineGC,
				$lineLetters,
				$lineBases;
			printf "%.3f\t", $lineGC / $lineLetters;
			printf "%.3f\n", $lineGC / $lineBases;
		}
	}
}
# update seq stats (seq name, GC content, seq length)
push @seqNamesList, substr($thisseq, 1);
push @seqLettersList, $seqLetters;
push @seqBasesList, $seqBases;
push @seqGCContentList, $seqGC;

print "\n" if $printConcatenatedGCStats;

$Nbases = $ACGTN[0] + $ACGTN[1] + $ACGTN[2] + $ACGTN[3];
$Nletters = $Nbases + $ACGTN[4];

# general GC content stats

if ($printTotalGCStats) {

	print "\n============================================================\n" if $verbose;
    print "Total GC content stats\n\n" if $verbose;

    print "statistic\tvalue\n";
    print "N seqs\t$Nseqs\n";
    print "N FASTA lines\t$Nfastalines\n";
    print "N lines in longest sequence\t$maxseqline\n";
    printf "Mean block length including N\t%.2f\n", $Nletters / $Nfastalines;
    printf "Mean block length excluding N\t%.2f\n", $Nbases / $Nfastalines;
    print "N letters including N\t$Nletters\n";
    print "N bases excluding N\t$Nbases\n";
    printf "GC including N\t%.5f\n", ($ACGTN[1] + $ACGTN[2]) / $Nletters;
    printf "GC excluding N\t%.5f\n", ($ACGTN[1] + $ACGTN[2]) / $Nbases;
    printf "A\t%d\n", $ACGTN[0];
    printf "C\t%d\n", $ACGTN[1];
    printf "G\t%d\n", $ACGTN[2];
    printf "T\t%d\n", $ACGTN[3];
    printf "N\t%d\n", $ACGTN[4];
    printf "A%%\t%.5f\n", $ACGTN[0] / $Nletters;
    printf "C%%\t%.5f\n", $ACGTN[1] / $Nletters;
    printf "G%%\t%.5f\n", $ACGTN[2] / $Nletters;
    printf "T%%\t%.5f\n", $ACGTN[3] / $Nletters;
    printf "N%%\t%.5f\n", $ACGTN[4] / $Nletters;

    print "\n";

}


# means across blocks in all sequences

if ($printBlockGCStats) {

	print "\n============================================================\n" if $verbose;
    print "means across blocks in all seqs\n\n" if $verbose;

    #print "seqGC: ", join(":", @seqGC), "\n";
    #print "seqLetters: ", join(":", @seqGC), "\n";
    #print "seqBases: ", join(":", @seqGC), "\n";
    print "block\tGC\tletters\tbases\tGCinclN\tGCexclN\n";
    foreach $i (0 .. $#seqGC) {
        printf "%d\t%d\t%d\t%d\t", 
                    $i, 
            $seqGC[$i],
            $seqLetters[$i],
            $seqBases[$i];
        printf "%.5f\t", $seqGC[$i] / $seqLetters[$i];
        printf "%.5f\n", $seqGC[$i] / $seqBases[$i];
    }
    print "\n";

}


# GC content for each sequence

if ($printSeqGCStats) {

	print "\n============================================================\n" if $verbose;
    print "GC content for each sequence\n\n" if $verbose;

    print "seq\tGC\tNletters\tNbases\tGCinclN\tGCexclN\n";
    foreach $i (0 .. $#seqNamesList) {
        print "$seqNamesList[$i]\t";
        print "$seqGCContentList[$i]\t";
        print "$seqLettersList[$i]\t";
        print "$seqBasesList[$i]\t";
        if ($seqLettersList[$i] == 0) {
            print "NA"
        } else {
            printf "%.5f", $seqGCContentList[$i] / $seqLettersList[$i];
        }
        print "\t";
        if ($seqBasesList[$i] == 0) {
            print "NA"
        } else {
            printf "%.5f", $seqGCContentList[$i] / $seqBasesList[$i];
        }
        print "\n";
    }
}

