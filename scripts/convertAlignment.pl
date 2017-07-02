#!/usr/bin/env perl

########################
#
#  This script interconverts alignment formats.
#
########################

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw();

use Bio::SeqIO;
use Bio::AlignIO;

my $o_infile;
my $o_informat = 'clustalw';
my $o_outfile;
my $o_outformat = 'fasta';
my $o_degap = 0;
my $o_uppercase = 0;
my $o_verbose = 0;
my $o_help = 0;

sub usage {
    print STDERR "
Usage: $0 [ --if input-format ] [ [--input] input.file ] [ --of output-format ] [ --output output.file ] 

Sequence format conversion - format to format

   --input     infile    default from STDIN, or first non-option command-line argument
   --output    outfile   default to STDOUT

   --if|format fasta     input format, see Bio::AlignIO documentation for available formats
               clustalw  http://bioperl.org/howtos/AlignIO_and_SimpleAlign_HOWTO.html
               nexus
               phylip    default $o_informat
               etc.

   --of|ormat  fasta     output file format, as above for formats
               etc.      default $o_outformat

   --degap               convert aligned sequences to standard sequences by removing gaps
   --uc | --uppercase    convert sequence to uppercase


   -v | --verbose        print end summary of conversion effort
   -h | --help           this help

";
    exit(1);
}

GetOptions (
    "input=s"      => \$o_infile,
    "output=s"     => \$o_outfile,
    "if|format=s"  => \$o_informat,
    "oformat=s"    => \$o_outformat,
    "degap"        => \$o_degap,
    "uc|uppercase" => \$o_uppercase,
    "verbose"      => \$o_verbose,
    "help"         => \$o_help
) or usage();

usage() if $o_help;

die("Only one input file allowed") if scalar(@ARGV) > 1;
$o_infile = $ARGV[0] if $ARGV[0];
die("Must specify input format and output format") if not $o_informat or not $o_outformat;

my ($in, $out);

$in  = ($o_infile)  ? Bio::AlignIO->new( -file => "<$o_infile", -format => $o_informat )
                    : Bio::AlignIO->new( -fh => \*STDIN, -format => $o_informat );

$out = ($o_outfile) ? Bio::AlignIO->new( -file => ">$o_outfile", -format => $o_outformat )
                    : Bio::AlignIO->new( -fh => \*STDOUT, -format => $o_outformat );

my $i = 0;

while (my $aln = $in->next_aln()) {
    ++$i;
    $aln->set_displayname_flat;
    if ($o_degap) {
        my $gc = $aln->gap_char();
        $aln->map_chars($gc, "");
    }
    $aln->uppercase() if $o_uppercase;
    $out->write_aln($aln);
}

print STDERR "$0 total of $i alignments converted from $o_informat to $o_outformat\n" if $o_verbose;

$in->close;
$out->close;


