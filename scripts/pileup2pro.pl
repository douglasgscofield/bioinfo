#! /usr/bin/perl -w

# Copyright (c) 2012 Douglas G. Scofield, Umeå Plant Sciences Centre, Umeå, Sweden
# douglas.scofield@plantphys.umu.se
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
use POSIX qw/isdigit/;
use Getopt::Long;

my $in_file;
my $out_file;
# contents of some pileup columns
my $seq_column = 0;
my $coord_column = 1;
my $ref_base_column = 2;
my $cvrg_column = 3;
my $read_bases_column = 4;
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

  samtools mpileup -B -C50 -q1 -f ref.fa your.bam | pileup2pro.pl > mlRho-input.txt

OPTIONS

    -                        read input from STDIN, write to STDOUT
    -in <filename>           read input from <filename>, else from STDIN
    -out <filename>          write output to <filename>, else to STDOUT
    -quiet                   don't print progress to STDERR
    -help, -?                help message

  Profile format lists bases present at each position in a reference sequence
  This script simply converts the format, so any filtering on base or mapping 
  quality, etc. that you may wish to do should be done when generating the pileup.

  mlRho estimates population genetic parameters from NGS data sequenced from a 
  diploid genome using Lynch's (2008) ML methods.  References are

  Haubold B, P Pfaffelhuber, and M Lynch. 2010. mlRho - a program for estimating
      the population mutation and recombination rates from shotgun-sequenced
      diploid genomes.  Molecular Ecology 19 Supplement 1:277-284.

  http://guanine.evolbio.mpg.de/mlRho

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
    "quiet" => \$quiet,
    "help|?" => \$help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $help;

if ($stdio or ! (defined($in_file) and ! $ARGV[0])) {
    *IN = *STDIN;
} else {
    $in_file = $ARGV[0] if ! $in_file;
    open (IN, "<$in_file") or die "Cannot open $in_file $!\n";
}

if ($stdio or ! defined($out_file)) {
    *OUT = *STDOUT;
} else {
    open (OUT, ">$out_file") or die "Cannot open $out_file $!\n";
}

while (<IN>) {
    #chomp;
    next if m/^\#/;
    my @fields = split /\t/;
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

    my $read_bases = $fields[ $read_bases_column ];
    die "Coverage column" . ($cvrg_column+1) . " contains non-numeric values. Check your input parameters as well as format of input dataset." if ( not isdigit $fields[ $cvrg_column ] );
    if ($read_bases =~ m/[\$\^\+-]/) {
        $read_bases =~ s/\^.//g; # removing the start of the read segment mark
        $read_bases =~ s/\$//g; # removing end of the read segment mark
        while ($read_bases =~ m/([\+-]){1}(\d+)/g) {
            my $indel_len = $2;
            $read_bases =~ s/([\+-]{1}$indel_len.{$indel_len})//; # remove indel info from read base field
        }
    }

    my @bases = split //, $read_bases;

    for my $base ( 0 .. @bases - 1 ) {
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

