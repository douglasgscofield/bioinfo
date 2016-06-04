#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use POSIX;

my $o_fai;
my $o_width = 10000;
my $o_overlap = 0;
my $o_nopad = 0;
my $o_noaddone = 0;
my $o_dropunderwidth = 0;

my $usage = "
SYNOPSIS

    $0 --fai file.fa.fai --width 100000 > file.windows100000.bed

Create BED file of evenly-spaced windows from Fasta index.

This script creates evenly-spaced windows along sequences described in a Fasta
index file ('.fai').  It is useful for creating BED files for sliding window
analyses in other tools.

Window width is set with --width, and window overlap (if any) is set with
--overlap.  The BED intervals are named after the sequence and the coordinate
of the beginning of the window plus 1 to match the convention in non-BED files,
i.e.,

      chr01 \\t 0 \\t 100000 \\t chr01_00000001
      chr01 \\t 100000 \\t 200000 \\t chr01_00100001

The last window for each sequence is shortened, if necessary, not to extend
beyond the end of the sequence; underwidth windows can be dropped from the
output by specifying the --drop-underwidth option.

The BED file of windows is written to stdout with sequences appearing in their
order within the '.fai' file and BED intervals in increasing order of starting
coordinate.

OPTIONS

    --fai FILE         REQUIRED: Fasta index file; this can be created with
                       'samtools faidx file.fa' and other tools
    --width INT        Window width [$o_width]
    --overlap INT      Window overlap [$o_overlap]
    --no-pad           Do not add 0-padding to the digits used to write the
                       starting coordinate in the interval name; without this,
                       the number of digits is calculated from the maximum
                       needed to make all coordinates equal-width in the
                       intervals across each sequence.
    --no-add-one       Do not add 1 to the interval name 
    --drop-underwidth  Do not include underwidth windows in the output; this
                       drops the last window for each sequence if it is not
                       the specified width.

";
sub dieusage { my $msg = join("\n", @_); die "\nERROR: $msg\n$usage"; }


GetOptions("fai=s"           => \$o_fai,
           "width=i"         => \$o_width,
           "overlap=i"       => \$o_overlap,
           "no-pad"          => \$o_nopad,
           "no-add-one"      => \$o_noaddone,
           "drop-underwidth" => \$o_dropunderwidth) or dieusage("unrecognized option");
dieusage("Must supply Fasta index file using --fai") if not $o_fai;
dieusage("--width must be > 0 and --overlap must be >= 0") if $o_width <= 0 or $o_overlap < 0;

my $pad = 0;
open (my $fai, "<", $o_fai) or die "could not open $o_fai: $!";
for (<$fai>) {
    chomp;
    my ($seqname, $seqlen, undef, undef, undef) = split /\t/;
    $pad = ceil(log10($seqlen)) if not $o_nopad;
    for (my $pos = 0; $pos < $seqlen; $pos += ($o_width - $o_overlap)) {
        my $w_name = "${seqname}_".sprintf("%0${pad}d", $pos + !$o_noaddone);
        my $w_end = $pos + $o_width;
        if ($w_end > $seqlen) {
            next if $o_dropunderwidth;
            $w_end = $seqlen;
        }
        print "$seqname\t$pos\t$w_end\t$w_name\n";
    }
}

