#!/usr/bin/env perl

########################
#
#  This script takes a sequence file and converts it to another sequence format
#
########################

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

my $o_infile = '';
my $o_informat = 'genbank';
my $o_outfile = '';
my $o_outformat = 'fasta';
my $o_verbose = 0;
my $o_help = 0;

sub usage {
    print STDERR "
Usage: $0 [ --if input-format ] [ [--input] input.file ] [ --of output-format ] [ --output output.file ] 

Sequence format conversion - format to format

   --input     infile    default from STDIN, or first non-option command line argument
   --output    outfile   default to STDOUT

   --if|format fasta     input format, see Bio::SeqIO documentation for available formats
               genbank   http://bioperl.org/howtos/SeqIO_HOWTO.html
               swiss
               embl      default $o_informat
               etc.

   --of|ormat  fasta     output file format, as above for formats
               etc.      default $o_outformat

   -v|--verbose          print end summary of conversion effort
   -h|--help             this help

";
    exit(1);
}

GetOptions (
    "input=s" => \$o_infile,
    "output=s" => \$o_outfile,
    "if|format=s" => \$o_informat,
    "oformat=s" => \$o_outformat,
    "verbose:1" => \$o_verbose,
    "help" => \$o_help,
) or usage();

usage() if $o_help;

die("Only one input file allowed") if scalar(@ARGV) > 1;
$o_infile = $ARGV[0] if $ARGV[0];
die("Must specify input format and output format") if not $o_informat or not $o_outformat;

my $in  = ($o_infile)  ? Bio::SeqIO->new( -file => "<$o_infile", -format => $o_informat )
                       : Bio::SeqIO->new( -fh => \*STDIN, -format => $o_informat );

my $out = ($o_outfile) ? Bio::SeqIO->new( -file => ">$o_outfile", -format => $o_outformat )
                       : Bio::SeqIO->new( -fh => \*STDOUT, -format => $o_outformat );

my $i = 0;

while (my $seq = $in->next_seq()) {
    ++$i;
    $out->write_seq($seq);
}

print STDERR "$0 total of $i sequences converted from $o_informat to $o_outformat\n" if $o_verbose;

$in->close;
$out->close;


