#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $names_file = "";
my $in_file = "";
my $out_file = "";
my $flag_contigmatched = 0;
my $Nmatched = 0;
my $Nnames = 0;
my $Nseqs = 0;
my $Nlines = 0;
my $Nletters = 0;
my $Nletters_matched = 0;
my $stdio;
my $help;

my $usage = "
NAME

    extractFastaSeqs.pl - extract named FASTA sequences

SYNOPSIS

    extractFastaSeqs.pl --names subset-names.txt --in full.fa --out subset.fa
    extractFastaSeqs.pl -n subset-names.txt -i full.fa - > subset.fa
    cat full.fa | extractFastaSeqs.pl -n subset-names.txt - > subset.fa
    extractFastaSeqs.pl subset-names.txt full.fa subset.fa

    Names must be given one per line in the names file.  Names of sequences 
    in the FASTA file are any characters after the initial '>' character in 
    the header of a FASTA sequence, so FASTA sequence names in the header as
    interpreted here are only delimited by end-of-line.

OPTIONS

    -                       read FASTA sequences from stdin and/or write
                            extracted sequences to stdout
    -n FILE, --names FILE   file containing names of FASTA sequences to extract
    -i FILE, --in FILE      input FASTA sequences
    -o FILE, --out FILE     output FASTA sequences

    -?, --help              help message

";

sub print_usage_and_exit($) {
    my $code = shift;
    print $usage;
    exit $code;
}

print_usage_and_exit(1) if ! defined($ARGV[0]);

GetOptions(
    "" => \$stdio, 
    "in=s" => \$in_file,
    "out=s" => \$out_file,
    "names=s" => \$names_file,
    "help|?" => \$help,
) or print_usage_and_exit(1);

print_usage_and_exit(0) if $help;

if (! $names_file and ! $in_file and ! $out_file) {
    $names_file = $ARGV[0];
    $in_file = $ARGV[1];
    $out_file = $ARGV[2];
}

$in_file = "/dev/stdin" if ! $in_file and $stdio;
$out_file = "/dev/stdout" if ! $out_file and $stdio;

print_usage_and_exit(1) if ! $names_file or ! $in_file or ! $out_file;

# print "names_file = $names_file, input = $in_file, output = $out_file\n";

open(NAMES, "<$names_file") or die("couldn't open names file: $!");
open(FASTA, "<$in_file") or die("couldn't open FASTA input: $!");
open(OUT, ">$out_file") or die("couldn't open FASTA output: $!");

my %NAMES;

while (<NAMES>) {
    ++$Nnames;
    chomp;
    print STDERR "$_ seen again in names file\n" if defined $NAMES{$_};
    ++$NAMES{$_};
}

while (<FASTA>) {
    ++$Nlines;
    chomp;
    my $linelength = 0;
    if (/^>/) {
        ++$Nseqs;
        if (defined($NAMES{substr($_, 1)})) {
            $flag_contigmatched = 1;
            ++$Nmatched;
        } else {
            $flag_contigmatched = 0;
        }
    } else {
        $linelength = length;
        $Nletters += $linelength;
    }
    if ($flag_contigmatched) {
        $Nletters_matched += $linelength;
        print OUT "$_\n";
    }
}
print STDERR "$Nseqs FASTA sequences in $Nlines lines, total length $Nletters bp\n";
print STDERR "$Nnames names, of these $Nmatched matched FASTA sequences, total length $Nletters_matched bp\n";

