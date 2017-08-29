#!/usr/bin/env perl

# Copyright (c) 2012,2015 Douglas G. Scofield, Uppsala University
# douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.

use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use File::Basename qw();
my $script = File::Basename::basename($0);

my $o_namesfile = "";
my $o_infile = "";
my $o_outfile = "";
my $o_header = 0;
my $o_inverse = 0;
my $o_quitonseen = 0;
my $o_stdio;
my $o_help;
my @o_names;

my $contigmatched = 0;
my $Nmatched = 0;
my $Nnames = 0;
my $Nseqs = 0;
my $Nlines = 0;
my $Nletters = 0;
my $Nletters_output = 0;

my $usage = "
NAME

    extractFastaSeqs.pl - extract named FASTA sequences

SYNOPSIS

    extractFastaSeqs.pl -n seq1 -n seq2 --in full.fa --out subset.fa
    extractFastaSeqs.pl --namesfile subset-names.txt --in full.fa --out subset.fa
    extractFastaSeqs.pl -nf subset-names.txt -i full.fa - > subset.fa
    cat full.fa | extractFastaSeqs.pl -nf subset-names.txt - > subset.fa
    extractFastaSeqs.pl subset-names.txt full.fa subset.fa

    For the -nf/--namesfile options, names must be given one per line in the
    names file.  Names of sequences in the FASTA file are any characters after the
    initial '>' character in the header of a FASTA sequence, followed by whitespace
    or an end of line, and must be matched in their entirety.  Use -H/--header to
    match against the entire header of the FASTA sequence.  Use -v/--inverse to
    extract sequences that _do not match_ any of the names.

    The --name and --namesfile options may be combined.

OPTIONS

    -                        read FASTA sequences from stdin and/or write
                             extracted sequences to stdout
    -n NAME, --name NAME     name of FASTA sequences to extract; may be specified
                             multiple times
    -nf FILE, --namesfile FILE
                             FILE contains names of FASTA sequences to extract
    -i FILE, --in FILE       input FASTA sequences
    -o FILE, --out FILE      output FASTA sequences
    -H, --header             match entire contents of the FASTA header
    -v, --inverse            output FASTA sequences that _do not match_ any of
                             names given.
    --quit-on-seen           quit once all sequences in the names file are seen

    -?, --help               help message
";

sub print_usage_and_exit($) {
    my $code = shift;
    say $usage;
    exit $code;
}

print_usage_and_exit(1) if ! defined($ARGV[0]);

GetOptions(
    "" => \$o_stdio, 
    "i|in=s" => \$o_infile,
    "o|out=s" => \$o_outfile,
    "nf|namesfile=s" => \$o_namesfile,
    "n|name=s" => \@o_names,
    "H|header" => \$o_header,
    "v|inverse" => \$o_inverse,
    "quit-on-seen" => \$o_quitonseen,
    "help|?" => \$o_help,
) or print_usage_and_exit(1);

print_usage_and_exit(0) if $o_help;

if (! $o_namesfile and ! $o_infile and ! $o_outfile) {
    $o_namesfile = $ARGV[0];
    $o_infile = $ARGV[1];
    $o_outfile = $ARGV[2];
}

$o_infile = "/dev/stdin" if ! $o_infile and $o_stdio;
$o_outfile = "/dev/stdout" if ! $o_outfile and $o_stdio;

print_usage_and_exit(1) if (! @o_names and ! $o_namesfile) or ! $o_infile or ! $o_outfile;

my %NAMES;

sub add_NAME($) {
    my $n = shift;
    if (defined $NAMES{$n}) {
        say STDERR "$script: sequence name $n specified again";
    } else {
        ++$Nnames;
        $NAMES{$n} = 0;
    }
}

add_NAME($_) foreach @o_names;

if ($o_namesfile) {
    open(my $NAMES, "<", $o_namesfile) or die("couldn't open names file: $!");
    while (<$NAMES>) {
        chomp;
        add_NAME($_);
    }
}

open(my $FASTA, "<", $o_infile) or die("couldn't open FASTA input: $!");
open(my $OUT, ">", $o_outfile) or die("couldn't open FASTA output: $!");

my $Nnames_seen = 0;

while (<$FASTA>) {
    ++$Nlines;
    chomp;
    my $linelength = 0;
    if (/^>/) {
        my $header = substr($_, 1);
        my @header_fields = split /\s\s*/, $header, 2;
        my $name_to_match = $o_header ? $header : $header_fields[0];
        if (defined($NAMES{$name_to_match})) {
            $contigmatched = 1;
            ++$Nmatched;
            ++$Nnames_seen if ($NAMES{$name_to_match} == 0);
            ++$NAMES{$name_to_match};
        } else {
            $contigmatched = 0;
            if ($o_quitonseen and $Nnames_seen == $Nnames) {
                --$Nlines;  # take back this line
                last;
            }
        }
        ++$Nseqs;
    } else {
        $linelength = length;
        $Nletters += $linelength;
    }
    if ($o_inverse xor $contigmatched) {
        $Nletters_output += $linelength;
        say $OUT "$_";
    }
}
say STDERR "$Nseqs FASTA sequences in $Nlines lines, total length $Nletters bp";
say STDERR "$Nnames unique names, and these matched $Nmatched FASTA sequences";
say STDERR "Output $Nletters_output bp";
say STDERR "Quit after seeing all $Nnames_seen names" if $o_quitonseen;

