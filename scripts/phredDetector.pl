#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);

my $opt_reads = 10000;
my $opt_verbose = 0;
my $opt_stdin = 0;
my $opt_help = 0;
my $infile = "";
my $opt_wide = 50;
my $opt_narrow = 20;

sub usage () {
   print STDERR "
Usage: 

    $0  [ options ]  fastq_file1.fq[.gz]

Input must be in FastQ format, if a filename is given it may be gzipped (*.gz)

Output to stdout is either '33', '64' or '59', depending on whether
Phred-scaled quality scores of base 33, base 64, or Solexa base 64 (beginning
with 59 ';') are detected.  The only data interpreted are read quality scores
on line 4 of each read; all sequence, pairing information etc. is ignored.

If all bases on input are of unusually high quality, then a Phred base of 59 or
64 may be reported when a Phred scale based on 33 was the one actually used.  A
few heuristics are used to detect possible problems, but these are not
comprehensive.

* If (maximum quality - minimum quality) >= $opt_wide, a warning message is printed
  to stderr and detected quality encoding (possibly erroneous) to stdout

* If (maximum quality - minimum quality) <= $opt_narrow, a warning message is printed
  to stderr and detected quality encoding (possibly erroneous) to stdout

* If otherwise unusual quality scores or unknown input were detected, an error
  message is printed to stderr and '??' to stdout

This script does not diagnose faulty FastQ files, nor does it fully diagnose
the various versions of Solexa, Sanger, Illumina pipelines described in
http://en.wikipedia.org/wiki/FASTQ_format.  It simply applies a minimum cutoff
of 33, 59 or 64 for quality values.  It is thus compatible with Sanger encoding but
does not take into account finer distinctions such as Illumina 1.3+ pipelines
not typically producing quality scores 0 and 1.

Options:

    -             Read uncompressed FastQ from stdin
    --reads INT   Number of reads to process to determine Phred basis [$opt_reads]
                  If 0, process *all* reads in the input file
    --wide INT    Use INT for the 'too wide' first heuristic above [$opt_wide]
    --narrow INT  Use INT for the 'too narrow' second heuristic above [$opt_narrow]

    --help | -?   Generate this help output

";
    exit 1;
}

GetOptions(
    "" => \$opt_stdin, 
    "reads" => \$opt_reads, 
    "wide" => \$opt_wide, 
    "narrow" => \$opt_narrow, 
    "help|?" => \$opt_help) or usage();
usage() if $opt_help or (!$ARGV[0] and !$opt_stdin);

if ($opt_stdin) {
    *INFILE = *STDIN;
} else {
    $infile = $ARGV[0];
    open(INFILE, "gzip -f -c -d ${infile} |") or die "couldn't open input $infile: $!";
}

my $min_quality = 1000;  # minimum quality seen so far
my $max_quality = -1000;  # maximum quality seen so far
my $n = 1;

while (<INFILE> and (!$opt_reads or $n <= $opt_reads)) {
    # Read name, beginning with '@'
	$_ = <INFILE>; # Sequence
	$_ = <INFILE>; # Optionally read name again, line begins with '+' no matter what
    $_ = <INFILE>; # Quality
    chomp;
    my @q = unpack('U*', $_);
    my $min = min(@q);
    my $max = max(@q);
    $min_quality = $min if $min < $min_quality;
    $max_quality = $max if $max > $max_quality;
    ++$n;
}
close(INFILE);

# Based on the above-cited wiki page
#
# if min_quality >= 33 ('!') then straight Sanger or Illumina 1.3 to pre-1.5
#
# if min_quality >= 59 (';') then Solexa
#
# if min_quality >= 66 ('B') then Illumina 1.5 to pre-1.8, note that 64 and 65 '@' 'A'
# are not produced by these pipelines
#
# if min_quality >= 35 ('#') then Illumina 1.8+, note that 33 and 34 '!' '"' are not 
# produced by these pipelines
#

my $range = $max_quality - $min_quality;
print STDERR "Warning: Phred quality range $range ('".chr($min_quality)."' to '".chr($max_quality)."') <= $opt_narrow, seems too narrow\n" if $range <= $opt_narrow;
print STDERR "Warning: Phred quality range $range ('".chr($min_quality)."' to '".chr($max_quality)."') >= $opt_wide, seems too wide\n" if $range >= $opt_wide;
if ($min_quality < 33 or $max_quality > 126) {
    print STDOUT "??\n";
} elsif ($min_quality >= 64) {
    print STDOUT "64\n";
} elsif ($min_quality >= 59) {
    print STDOUT "59\n";
} elsif ($min_quality >= 33) {
    print STDOUT "33\n";
}

