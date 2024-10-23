#!/usr/bin/env perl

# Produce BED file(s) describing reference sequences from a .fai-format Fasta index
#
# Copyright (c) 2015 Douglas G. Scofield, Uppsala University
# douglasgscofield@gmail.com
# 
# No warranty is implied or should be assumed with this code.
#
# CHANGELOG
#
# 2015-10-13
# - create the script
#
# TODO

use strict;
use warnings;
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
my $list_file;
my $N_list;
my $stdio;
my $o_quiet = 0;
my $help;
my $o_header = 0;
my $N_references = 0; # the number of reference sequences seen

my $usage = "
NAME

  fai2Bed.pl - create BED file from .fai-format Fasta index

SYNOPSIS

  fai2Bed.pl --min-length 1000 fasta.fa.fai > min1000.bed

  fai2Bed.pl --list seqs.txt fasta.fa.fai > seqs-from-fasta.bed

EXAMPLE

OPTIONS

    -                     Read .fai file from stdin, write to stdout
    -i FILE, --in FILE    Read .fai file from FILE, else from stdin
    -o FILE, --out FILE   Write BED output to FILE, else to stdout

    --num INT             Include just the first INT sequences
    --chunk-size INT      Produce separate BED files each containing sequences
                          representing approximately INT bp; reference sequences
                          are only described complete, so each BED is likely to 
                          describe more than INT bp.  Individual BED files are
                          named FILE.xx.bed where FILE is specified with --out,
                          which is required, and x is the integer sequence of 
                          file creation.
    --min-length INT      Do not include reference sequences shorter than INT
    --max-length INT      Do not include reference sequences longer than INT

    -L FILE, --list FILE  Include the sequences listed in FILE in the BED and no
                          others. Cannot be used with other selection options.

    --header              Include a BED header line in BED output

    -q, --quiet           Do not produce informational messages to stderr

    -?, --help            Help message

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
    "chunk-size|chunksize=f" => \$opt_chunksize,
    "min-length|minlength=f" => \$min_length,
    "max-length|maxlength=f" => \$max_length,
    "L|list=s" => \$list_file,
    "header" => \$o_header,
    "q|quiet" => \$o_quiet,
    "help|?" => \$help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $help;
print_usage_and_exit("") if $opt_chunksize < 0;
print_usage_and_exit("") if $opt_chunksize and not $out_file;
print_usage_and_exit("") if $min_length > $max_length;
print_usage_and_exit("") if $list_file and ($min_length || $max_length != 9.99e12 || $num || $opt_chunksize);

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
print $OUT "#chrom\tchromStart\tchromEnd\n" if $o_header;

my $LIST = FileHandle->new;
my %LIST;
if ($list_file) {
    $LIST->open("< $list_file") or die "Cannot open $list_file: $!\n";
    while (<$LIST>) {
        next if /^#/;
        chomp;
        $LIST{$_}++;
    }
    undef $LIST;
    $N_list = scalar keys %LIST;
    die "no sequences in list file '$list_file'" if $N_list == 0;
}

while (<$IN>) {
    # Format of an .fai line is (https://www.biostars.org/p/11523/#11524)
    # 1. sequence name
    # 2. sequence length
    # 3. offset of the first base of the chromosome sequence in the file
    # 4. number of bases in each fasta line
    # 5. number of bytes in each fasta line
    chomp;
    my @fields = split /\t/;
    my $ref = $fields[0];
    my $length = $fields[1];

    if ($list_file && ! defined $LIST{$ref}) {
        next;
    } elsif ($length < $min_length or $length > $max_length) {
        next;
    }

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
        print $OUT "#chrom\tchromStart\tchromEnd\n" if $o_header;
    }

    print $OUT "$ref\t0\t$length\n";

    if ($opt_chunksize and $this_chunk_size > $opt_chunksize) {
        print STDERR "$out_file represents $this_chunk_seqs sequences containing $this_chunk_size bp\n" if ! $o_quiet;
        $OUT->close();
        $this_chunk_size = $this_chunk_seqs = 0;
        $must_open_chunk = 1;  # open a new bed file
    }
}

print STDERR "$out_file represents $this_chunk_seqs sequences containing $this_chunk_size bp\n" if ! $o_quiet && ! $must_open_chunk;
print STDERR "$list_file contained $N_list sequence names\n" if $list_file;

undef $OUT;
undef $IN;


