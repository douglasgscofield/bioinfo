#! /usr/bin/perl -w

# Produce BED file(s) describing reference sequences from a SAM header
#
# Copyright (c) 2012 Douglas G. Scofield, Umeå Plant Sciences Centre, Umeå, Sweden
# douglasgscofield@gmail.com
# 
# No warranty is implied or should be assumed with this code.
#
# CHANGELOG
#
# 2012-11-20
# - print out stats for final BED, fixups to option variables and 
#   GetOptions call 
# 2012-11-10
# - add --num and --chunk-size options, change option names, switch to
#   new style file handles, turning the code into a bit of a mess due to
#   need for rapid implementation
# 2012-05-19
# - create the script
#
# TODO

use strict;
use FileHandle;
use Getopt::Long;

my $min_length = 0;
my $max_length = 9.99e12;
my $num = 0;
my $opt_chunksize = 0;
my $N_chunk = 0;
my $chunk_prefix;
my $this_chunk_size = 0; # the total length of this chunk
my $this_chunk_seqs = 0; # the total number of sequences in this chunk

my $in_file;
my $out_file;
my $stdio;
my $quiet = 0;
my $help;
my $opt_header = 0;
my $N_references = 0; # the number of reference sequences seen

my $usage = "
NAME

  samHeader2bed.pl - create BED file from SAM header containing length-filtered contigs

SYNOPSIS

  samtools view -H your.bam | samHeader2bed.pl -min-length 1000 - > min1000.bed
  samtools mpileup -l min1000.bed -u -f ref.fa your.bam | bcftools view ...

EXAMPLE

OPTIONS

    -                     read SAM header from stdin, write to stdout
    -i FILE, --in FILE    read SAM header from FILE, else from stdin
    -o FILE, --out FILE   write BED output to FILE, else to stdout
    --num INT             include just the first INT sequences
    --chunk-size INT      produce separate BED files each containing sequences
                          representing approximately INT bp; reference sequences
                          are only described complete, so each BED is likely to 
                          describe more than INT bp.  Individual BED files are
                          named FILE.xx.bed where FILE is specified with --out,
                          which is required, and x is the integer sequence of 
                          file creation
    --min-length INT      do not include reference sequences shorter than INT
    --max-length INT      do not include reference sequences longer than INT
    --header              include a BED header line in BED output

    -?, --help            help message

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
    "num=f" => \$num,
    "chunk-size=f" => \$opt_chunksize,
    "min-length|minlength=f" => \$min_length,
    "max-length|maxlength=f" => \$max_length,
    "header" => \$opt_header,
    "help|?" => \$help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $help;
print_usage_and_exit("") if $opt_chunksize < 0;
print_usage_and_exit("") if $opt_chunksize and not $out_file;
print_usage_and_exit("") if $min_length > $max_length;

my $IN = FileHandle->new;
$in_file = $ARGV[0] if $ARGV[0] and ! $in_file;
$in_file = "/dev/stdin" if ! defined($in_file) or $stdio;
$IN->open("< $in_file") or die "Cannot open $in_file: $!\n";

my $OUT = FileHandle->new;
my $must_open_chunk = 0;
if ($opt_chunksize) {
    $chunk_prefix = $out_file;
    $N_chunk = 1;
    $this_chunk_size = $this_chunk_seqs = 0;
    $out_file = sprintf("%s.%.2d.bed", $chunk_prefix, $N_chunk);
} else {
    $out_file = "/dev/stdout" if $stdio or ! defined($out_file);
}
$OUT->open("> $out_file") or die "Cannot open $out_file: $!\n";
print $OUT "#chrom\tchromStart\tchromEnd\n" if $opt_header;

while (<$IN>) {
    next if ! m/^\@SQ/;
    chomp;
    my @fields = split /\t/;
    my $ref = substr($fields[1], 3);  # everything after SN: in the field
    my $length = substr($fields[2], 3);  # everything after LN: in the field

    next if $length < $min_length or $length > $max_length;

    ++$N_references;

    last if $num and $N_references > $num;

    ++$this_chunk_seqs;
    $this_chunk_size += $length;

    if ($must_open_chunk) { # we need to open a new bed file for the next chunk
        $must_open_chunk = 0;
        ++$N_chunk;
        #$out_file = $chunk_prefix . ".$N_chunk";
        $out_file = sprintf("%s.%.2d.bed", $chunk_prefix, $N_chunk);
        $OUT->open("> $out_file") or die "Cannot open next chunk $out_file: $!\n";
        print $OUT "#chrom\tchromStart\tchromEnd\n" if $opt_header;
    }

    print $OUT "$ref\t0\t$length\n";

    if ($opt_chunksize and $this_chunk_size > $opt_chunksize) {
        print STDERR "chunk $out_file represents $this_chunk_seqs sequences containing $this_chunk_size bp\n";
        $OUT->close();
        $this_chunk_size = $this_chunk_seqs = 0;
        $must_open_chunk = 1;  # open a new bed file
    }
}

print STDERR "chunk $out_file represents $this_chunk_seqs sequences containing $this_chunk_size bp\n" if ! $must_open_chunk;

undef $OUT;
undef $IN;


