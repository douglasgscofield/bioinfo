#!/usr/bin/env perl

########################
#
#  This script takes a sequence file and converts it to another sequence format
#
########################

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(first);
use Bio::SeqIO;

sub usage {
    print STDERR <<__usage__;

Usage: $THISSCRIPT [ - ] [ input.file ] -f input-format [ -out output.file ] -of output-format 

Sequence format conversion - format to format

   infile
   -out|utput outfile

   -                    read input from stdin and/or write output to stdout

   -f|ormat fasta       input file format [default $informat], see
            genbank     Bio::SeqIO documentation for available formats
            swiss
            embl
            etc.

   -of  fasta           output file format, as above for input formats
        etc.

   -v|erbose            print end summary of conversion effort

__usage__
    exit(1);
}

GetOptions (
    "input=s" => \$infile,
    "output=s" => \$outfile,
    "format=s" => \$informat,
    "of=s" => \$outformat,
    "verbose:1" => \$verbose,
    "" => \$stdio
) or usage();

$infile = $ARGV[0] if $ARGV[0];
die("Only one input file allowed") if scalar(@ARGV) > 1;
die("Must specify input and/or output file without -") if !$stdio and !($infile and $outfile);
die("Must specify input format") if !$informat;
die("Must specify output format") if !$outformat;

my $in = ($infile) ? Bio::SeqIO->new( -file => "<$infile", -format => $informat)
                   : Bio::SeqIO->new( -fh => \*STDIN, -format => $informat);
my $out = ($outfile) ? Bio::SeqIO->new( -file => ">$outfile", -format => $outformat)
                     : Bio::SeqIO->new( -fh => \*STDOUT, -format => $outformat);

my $i = 0;

while (my $seq = $in->next_seq()) {
    ++$i;
    $out->write_seq($seq);
}

print STDERR "$THISSCRIPT: total of $i sequences converted from $informat to $outformat\n" if $verbose;

$in->close;
$out->close;


