#!/usr/bin/env perl

# Copyright (c) 2012,2015 Douglas G. Scofield, Uppsala University
# douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.
#
# Pileup-parsing code built upon the bones of the mpileup-format parser script
# https://bitbucket.org/galaxy/galaxy-central/src/tip/tools/samtools/pileup_parser.pl
#
# pileup2pro2.pl - Convert pileup to a profile format
#
# Profile2 is a format in which bases present at each position in a reference
# sequence are enumerated.  This script simply converts the format, so any
# filtering on base/mapping quality, etc. that you may wish to do should be
# done when generating the pileup:
#
#     samtools mpileup -B -C50 -q1 -f test.fa test.bam | pileup2pro2.pl > test_pro2.txt
#
# Profile2 format contains reference names, coordinates and raw numbers of bases
# (A, C, G, T, N), separated by tabs:
#
#     contig_1	1	0	2	0	0	0
#     contig_1	2	2	0	0	0	0
#     contig_1	3	0	2	0	0	0
#     contig_1	4	2	0	0	0	0
#     contig_1	5	0	0	0	2	0
#     contig_1	6	0	0	2	0	0
#     ...
#
# CHANGELOG
#
# 2015-05-20
# - create the script
#
# TODO

use strict;
use warnings;
use POSIX qw/isdigit/;
use Getopt::Long;

my $in_file;
my $out_file;
# contents of some pileup columns
my $seq_column = 0;
my $coord_column = 1;
my $ref_base_column = 2;
my $BAM_start_column = 3;

my $offset_coverage = 0;
my $offset_base_call = 1;
my $offset_base_quality = 2;

my @coverage_columns;
my @base_call_columns;
my @base_quality_columns;

my $which_bams = "";
my @which_bams = ();
my $o_stdio;
my $o_quiet;
my $o_help;
my $o_N;
my $current_reference = ""; # the name of the current reference sequence
my $N_references = 0; # the number of reference sequences seen
my $N_coordinates = 0; # the total number of coordinates seen

my $usage = "
NAME

  $0 - Convert pileup to profile2 format

SYNOPSIS

  samtools mpileup -B -C50 -q1 -f test.fa test.bam | pileup2pro2.pl > test_pro2.txt

OPTIONS

    -                          read input from stdin, write to stdout
    --in FILE                  read input from FILE, else from first non-argument
                               on the command line or from stdin
    --out FILE                 write output to FILE, else to stdout
    --which-bams INT[,INT...]  produce profile output for the INT-th BAM(s) in 
                               order as provided on the samtools mpileup command
                               line, starting with 1; otherwise produce profile
                               output for all BAMs 
    --N                        include 5th count column for Ns
    --quiet                    don't print progress to stderr
    --help, -?                 help message

  Profile2 format lists bases present at each position in a reference sequence,
  with columns sequence, position, A, C, G, T:

     contig_1	1	0	2	0	0
     contig_1	2	2	0	0	0
     contig_1	3	0	2	0	0
     contig_1	4	2	0	0	0
     contig_1	5	0	0	0	2
     contig_1	6	0	0	2	0

  This script simply converts the format so any filtering on base or mapping 
  quality, etc. that you may wish to do should be done when generating the pileup.
  Pileup format created from multiple BAM files has 3 columns per BAM file; this
  script will merge all columns while creating profile output up line unless the
  --which-bam option is given.

";

sub print_usage_and_exit($) {
    my $msg = shift;
    print "$msg\n" if $msg;
    print $usage;
    exit 0;
}

GetOptions(
    ""            => \$o_stdio, 
    "in=s"        => \$in_file,
    "out=s"       => \$out_file,
    "which-bam=s" => \$which_bams,
    "N"           => \$o_N,
    "quiet"       => \$o_quiet,
    "help|?"      => \$o_help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $o_help;

@which_bams = sort split(/,/, $which_bams) if $which_bams;

if ($o_stdio or ! (defined($in_file) or $ARGV[0])) {
    *IN = *STDIN;
} else {
    $in_file = $ARGV[0] if ! $in_file;
    open (IN, "<$in_file") or die "Cannot open $in_file: $!";
}

if ($o_stdio or ! defined($out_file)) {
    *OUT = *STDOUT;
} else {
    open (OUT, ">$out_file") or die "Cannot open $out_file: $!";
}

my $ncols = 0;

while (<IN>) {
    #chomp;
    next if m/^\#/;
    my @fields = split /\t/;
    if ($ncols and scalar(@fields) ne $ncols) {
    	print_usage_and_exit("inconsistent number of columns in pileup input line $.") if scalar(@fields) ne $ncols;
    }
    if (! $ncols) {
        # do setup after first line
        $ncols = scalar(@fields);
        my $n_bams = ($ncols - 3) / 3;
        @which_bams = 1 .. $n_bams if ! $which_bams;
        ((/\D/ or $_ < 1 or $_ > $n_bams) and die("must specify INT >= 1 and <= $n_bams")) foreach @which_bams;
        @which_bams = map { (($_ - 1) * 3) + 3 } @which_bams;
        @coverage_columns = map { $_ + $offset_coverage } @which_bams;
        @base_call_columns = map { $_ + $offset_base_call } @which_bams;
        @base_quality_columns = map { $_ + $offset_base_quality } @which_bams;
    }
    ++$N_coordinates;

    # do we have a new reference sequence?
    if ($current_reference ne $fields[ 0 ]) {
    	$current_reference = $fields[ 0 ];
        ++$N_references;
    }

    use constant SNPs_zeroes => ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'N' => 0);
    my %SNPs = SNPs_zeroes;
    if ($fields[ 2 ] eq "*") {
        print OUT $fields[0] . "\t" . $fields[1] . "\t-1\t-1\t-1\t-1".($o_N ? "\t-1" :"")."\n";
        #next; # skip indel lines
    }

    (/\D/ and print_usage_and_exit("error on input line $., coverage column invalid '$_'")) foreach @fields[ @coverage_columns ];
    my $read_bases = "";
    $read_bases .= $_ foreach @fields [ @base_call_columns ];  # concatenate base calls
    
    if ($read_bases =~ m/[\$\^*\+-]/) {
        $read_bases =~ s/\^.//g; # removing the start of the read segment mark
        $read_bases =~ s/[*\$]//g; # remove end of read mark and gap asterisks
        #$read_bases =~ s/\$//g; # removing end of the read segment mark
        #$read_bases =~ s/\*//g; # removing any asterisks that are there for gaps
        while ($read_bases =~ m/([\+-]){1}(\d+)/g) {
            my $indel_len = $2;
            $read_bases =~ s/([\+-]{1}$indel_len.{$indel_len})//; # remove indel info from read base field
        }
    }

    my @bases = split //, $read_bases;

    for my $base ( 0 .. $#bases ) {
        if ( $bases[ $base ] =~ m/[ATGCN]/i ) {
            $SNPs{ uc( $bases[ $base ] ) } += 1;
        } elsif ( $bases[ $base ] =~ m/[\.,]/ ) {
            $SNPs{ uc( $fields[ 2 ] ) } += 1;
        }           
    } 

    print OUT $fields[0] . "\t" . $fields[1] . "\t" .
              $SNPs{A} . "\t" . $SNPs{C} . "\t" . $SNPs{G} . "\t" .  $SNPs{T} . ($o_N ? ("\t".$SNPs{N}) : "") . "\n";

    print STDERR "seen $N_coordinates positions across $N_references reference sequences\n" if ! ($N_coordinates % 10000) and ! $o_quiet;
}

close IN;
close OUT;

