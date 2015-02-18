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

my $THISSCRIPT = File::Basename::basename($0);

my $stdio = 0;
my $infile = "";
my $outfile = "";
my $informat = "fasta";
my $outformat = "";
my $format = "";
my $verbose = 0;

sub usage {
    print STDERR <<__usage__;

Usage: $THISSCRIPT [ - ] [ input.file ] -f input-format [ -out output.file ] -of output-format 

Alignment conversion - format to format

   infile
   -out|utput outfile

   -                    read input from stdin and/or write output to stdout

   -f|ormat fasta       input file format [default $informat], see
            clustalw    Bio::AlignIO documentation for available formats
            nexus
            phylip
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

my $in = ($infile) ? Bio::AlignIO->new( -file => "<$infile", -format => $informat)
                   : Bio::AlignIO->new( -fh => \*STDIN, -format => $informat);
my $out = ($outfile) ? Bio::AlignIO->new( -file => ">$outfile", -format => $outformat)
                     : Bio::AlignIO->new( -fh => \*STDOUT, -format => $outformat);

my $i = 0;

while (my $aln = $in->next_aln()) {
    ++$i;
    $aln->set_displayname_flat;
    $out->write_aln($aln);
}

print STDERR "$THISSCRIPT: total of $i alignments converted from $informat to $outformat\n" if $verbose;

$in->close;
$out->close;


