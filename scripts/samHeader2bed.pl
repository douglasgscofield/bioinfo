#! /usr/bin/perl -w

# Produce a BED file containing length-filtered reference sequences from a SAM header
#
# Copyright (c) 2012 Douglas G. Scofield, Umeå Plant Sciences Centre, Umeå, Sweden
# douglasgscofield@gmail.com
# 
# No warranty is implied or should be assumed with this code.
#
# CHANGELOG
#
# 2012-05-19
# - create the script
#
# TODO

use strict;
use Getopt::Long;

my $minlength = 0;
my $maxlength = 9.99e12;

my $in_file;
my $out_file;
my $stdio;
my $help;
my $noheader = 0;
my $N_references = 0; # the number of reference sequences seen

my $usage = "
NAME

  samHeader2bed.pl - create BED file from SAM header containing length-filtered contigs

SYNOPSIS

  samtools view -H your.bam | samHeader2bed.pl -minlength 1000 - > min1000.bed
  samtools mpileup -l min1000.bed -u -f ref.fa your.bam | bcftools view ...

EXAMPLE

OPTIONS

    -                        read SAM header from STDIN, write to STDOUT
    -in <filename>           read SAM header from <filename>, else from STDIN
    -out <filename>          write BED output to <filename>, else to STDOUT
    -noheader                do not print header on output
    -minlength <length>      do not include reference sequences shorter than <length>
    -maxlength <length>      do not include reference sequences longer than <length>

    -help, -?                help message

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
    "noheader" => \$noheader,
    "minlength=f" => \$minlength,
    "maxlength=f" => \$maxlength,
    "help|?" => \$help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $help;

if ($stdio or (! defined($in_file) and ! $ARGV[0])) {
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

if (! $noheader) {
    print OUT "#chrom\tchromStart\tchromeEnd\n";
}

while (<IN>) {
    next if ! m/^\@SQ/;
    chomp;
    my @fields = split /\t/;
    ++$N_references;
    my $ref = substr($fields[1], 3);
    my $length = substr($fields[2], 3);

    if ($length >= $minlength and $length <= $maxlength) {
        print OUT "$ref\t0\t$length\n";
    }
}

close OUT;
