#!/usr/bin/env perl

# Copyright (c) 2012,2015 Douglas G. Scofield, Uppsala University
# douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.
#
# Pileup-parsing code built upon the bones of the mpileup-format parser script
# https://bitbucket.org/galaxy/galaxy-central/src/tip/tools/samtools/pileup_parser.pl
#
# pileup2pro.pl - Convert pileup to profile format as used for input to mlRho.  This
#       is a format in which bases present at each position in a reference sequence
#       are enumerated.  This script simply converts the format, so any filtering on 
#       base/mapping quality, etc. that you may wish to do should be done when 
#       generating the pileup:
#
#       samtools mpileup -B -C50 -q1 -f ref.fa your.bam | pileup2pro.pl > mlRho-input.txt
#
#       mlRho estimates population genetic parameters from NGS data sequenced from a 
#       diploid genome.
#
#       Haubold B, P Pfaffelhuber, and M Lynch. 2010. mlRho - a program for estimating
#           the population mutation and recombination rates from shotgun-sequenced
#           diploid genomes.  Molecular Ecology 19 Supplement 1:277-284.
#
#           http://guanine.evolbio.mpg.de/mlRho
#
#       It implements methods described in
#
#       Lynch M. 2008. Estimation of nucleotide diversity, disequilibrium coefficients, 
#           and mutation rates from high-coverage genome-sequencing projects. Molecular
#           Biology and Evolution 25:2421-2431. 
#
#       Profile format contains reference names, coordinates and raw numbers of bases:
#
#       >contig_1
#       1	0	2	0	0
#       2	2	0	0	0
#       3	0	2	0	0
#       4	2	0	0	0
#       5	0	0	0	2
#       6	0	0	2	0
#       ...
#
# CHANGELOG
#
# 2012-05-19
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
my $has_mapping_quality = 0;
my $stdio;
my $quiet;
my $help;
my $current_reference = ""; # the name of the current reference sequence
my $N_references = 0; # the number of reference sequences seen
my $N_coordinates = 0; # the total number of coordinates seen

my $usage = "
NAME

  pileup2pro.pl - Convert pileup to profile format as used for input to mlRho

SYNOPSIS

  samtools mpileup -B -q1 -f ref.fa your.bam | pileup2pro.pl > mlRho-input.txt

OPTIONS

    -                          read input from stdin, write to stdout
    --in FILE                  read input from FILE, else from stdin
    --out FILE                 write output to FILE, else to stdout
    --which-bams INT[,INT...]  produce profile output for the INT-th BAM(s) in 
                               order as provided on the samtools mpileup command
                               line, starting with 1; otherwise produce profile
                               output for all BAMs 
    --has-mapping-quality      must be specified if -s used for samtools mpileup
    --quiet                    don't print progress to stderr
    --help, -?                 help message

  Profile format lists bases present at each position in a reference sequence.
  This script simply converts the format so any filtering on base or mapping 
  quality, etc. that you may wish to do should be done when generating the pileup.
  Pileup format created from multiple BAM files has 3 columns per BAM file; this
  script will merge all columns while creating profile output up line unless the
  --which-bam option is given.

  mlRho (Haubold et al. 2010) estimates population genetic parameters from NGS 
  data for a diploid genome using Lynch's (2008) ML methods.

  http://guanine.evolbio.mpg.de/mlRho

  Haubold B, P Pfaffelhuber, and M Lynch. 2010. mlRho - a program for estimating
      the population mutation and recombination rates from shotgun-sequenced
      diploid genomes.  Molecular Ecology 19 Supplement 1:277-284.

  Lynch M. 2008. Estimation of nucleotide diversity, disequilibrium coefficients, 
      and mutation rates from high-coverage genome-sequencing projects. Molecular
      Biology and Evolution 25:2421-2431. 

";

sub print_usage_and_exit($) {
    my $msg = shift;
    print "$msg\n" if $msg;
    print $usage;
    exit 0;
}

GetOptions(
    "" => \$stdio, 
    "in=s" => \$in_file,
    "out=s" => \$out_file,
    "which-bam=s" => \$which_bams,
    "has-mapping-quality" => \$has_mapping_quality,
    "quiet" => \$quiet,
    "help|?" => \$help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $help;

@which_bams = sort split(/,/, $which_bams) if $which_bams;

if ($stdio or ! (defined($in_file) or $ARGV[0])) {
    *IN = *STDIN;
} else {
    $in_file = $ARGV[0] if ! $in_file;
    open (IN, "<$in_file") or die "Cannot open $in_file: $!";
}

if ($stdio or ! defined($out_file)) {
    *OUT = *STDOUT;
} else {
    open (OUT, ">$out_file") or die "Cannot open $out_file: $!";
}

my $ncols = 0;

sub field_mapping_quality_fixup(@)  # only call if $has_mapping_quality is given!!!
{
    my @fixup = @_[ $seq_column, $coord_column, $ref_base_column ];
    my $i = $BAM_start_column;
    while ($i <= $#_) {
	push @fixup, $_[ $i ];  # coverage
	push @fixup, $_[ $i + 1 ];  # base calls
	push @fixup, $_[ $i + 2 ];  # base quality
        # throw away mapping-quality column, not present with 0 coverage
        $i += (($_[ $i ] == 0) ? 3 : 4);
    }
    return(@fixup);
}

while (<IN>) {
    #chomp;
    next if m/^\#/;
    my @fields = split /\t/;
    if ($ncols and scalar(@fields) ne $ncols) {
        @fields = field_mapping_quality_fixup(@fields) if $has_mapping_quality;
    	print_usage_and_exit("inconsistent number of columns in pileup input line $., might you need --has-mapping-quality?") if scalar(@fields) ne $ncols;
    }
    if (! $ncols) {
        # do setup after first line
        @fields = field_mapping_quality_fixup(@fields) if $has_mapping_quality;
        $ncols = scalar(@fields);
        my $n_bams = ($ncols - 3) / 3;
        @which_bams = 1 .. $n_bams if ! $which_bams;
        ((/\D/ or $_ < 1 or $_ > $n_bams) and die("must specify INT >= 1 and <= $n_bams")) foreach @which_bams;
        @which_bams = map { (($_ - 1) * 3) + $BAM_start_column } @which_bams;
        @coverage_columns = map { $_ + $offset_coverage } @which_bams;
        @base_call_columns = map { $_ + $offset_base_call } @which_bams;
        @base_quality_columns = map { $_ + $offset_base_quality } @which_bams;
    }
    ++$N_coordinates;

    # do we have a new reference sequence?
    if ($current_reference ne $fields[ $seq_column ]) {
    	$current_reference = $fields[ $seq_column ];
        ++$N_references;
        print OUT ">$current_reference\n";
    }

    my @output_fields = ();
    use constant SNPs_zeroes => ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0);
    my %SNPs = SNPs_zeroes;
    if ($fields[ $ref_base_column ] eq "*") {
        next; # skip indel lines
    }

    (/\D/ and print_usage_and_exit("error on input line $., coverage column invalid '$_', might you need --has-mapping-quality?")) foreach @fields[ @coverage_columns ];
    my $read_bases = "";
    $read_bases .= $_ foreach @fields [ @base_call_columns ];  # concatenate base calls
    
    if ($read_bases =~ m/[\$\^*\+-]/) {
        $read_bases =~ s/\^.//g; # removing the start of the read segment mark
        $read_bases =~ s/\$//g; # removing end of the read segment mark
        $read_bases =~ s/\*//g; # removing any asterisks that are there for gaps
        while ($read_bases =~ m/([\+-]){1}(\d+)/g) {
            my $indel_len = $2;
            $read_bases =~ s/([\+-]{1}$indel_len.{$indel_len})//; # remove indel info from read base field
        }
    }

    my @bases = split //, $read_bases;

    for my $base ( 0 .. $#bases ) {
        if ( $bases[ $base ] =~ m/[ATGC]/i ) {
            $SNPs{ uc( $bases[ $base ] ) } += 1;
        } elsif ( $bases[ $base ] =~ m/[\.,]/ ) {
            $SNPs{ uc( $fields[ $ref_base_column ] ) } += 1;
        }           
    } 
    push @output_fields, $fields[ $coord_column ];

    foreach my $SNP (sort keys %SNPs) {
        push @output_fields, $SNPs{$SNP};
    }

    print OUT join("\t", @output_fields), "\n";

    print STDERR "seen $N_coordinates positions across $N_references reference sequences\n" if ! ($N_coordinates % 10000) and ! $quiet;
}

close IN;
close OUT;

