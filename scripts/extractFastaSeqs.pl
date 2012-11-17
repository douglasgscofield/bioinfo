#!/usr/bin/perl

# Copyright (c) 2012 Douglas G. Scofield, Umeå Plant Sciences Centre, Umeå, Sweden
# douglas.scofield@plantphys.umu.se
# douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw();
my $script = File::Basename::basename($0);

my $names_file = "";
my $in_file = "";
my $out_file = "";
my $opt_header = 0;
my $opt_reverse = 0;
my $opt_quitonseen = 0;
my $contigmatched = 0;
my $Nmatched = 0;
my $Nnames = 0;
my $Nseqs = 0;
my $Nlines = 0;
my $Nletters = 0;
my $Nletters_output = 0;
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
    the header of a FASTA sequence, followed by whitespace or an end of line,
    and must be matched in their entirety.  Use --header to match against the 
    entire header of the FASTA sequence.  Use --reverse to extract sequences 
    that _do not match_ any of the names.

OPTIONS

    -                       read FASTA sequences from stdin and/or write
                            extracted sequences to stdout
    -n FILE, --names FILE   file containing names of FASTA sequences to extract
    -i FILE, --in FILE      input FASTA sequences
    -o FILE, --out FILE     output FASTA sequences
    --header                match entire contents of the FASTA header
    --reverse               output FASTA sequences that _do not match_ any of
                            names given.
    --quit-on-seen          quit once all sequences in the names file are seen

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
    "header" => \$opt_header,
    "reverse" => \$opt_reverse,
    "quit-on-seen" => \$opt_quitonseen,
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
    chomp;
    if (defined $NAMES{$_}) {
        print STDERR "$script: $_ seen again in $names_file, line $.\n";
    } else {
        ++$Nnames;
        $NAMES{$_} = 0;
    }
}

my $Nnames_seen = 0;

while (<FASTA>) {
    ++$Nlines;
    chomp;
    my $linelength = 0;
    if (/^>/) {
        my $header = substr($_, 1);
        my @header_fields = split /\s\s*/, $header, 2;
        my $name_to_match = $opt_header ? $header : $header_fields[0];
        if (defined($NAMES{$name_to_match})) {
            $contigmatched = 1;
            ++$Nmatched;
            ++$Nnames_seen if ($NAMES{$name_to_match} == 0);
            ++$NAMES{$name_to_match};
        } else {
            $contigmatched = 0;
            if ($opt_quitonseen and $Nnames_seen == $Nnames) {
                --$Nlines;  # take back this line
                last;
            }
        }
        ++$Nseqs;
    } else {
        $linelength = length;
        $Nletters += $linelength;
    }
    if ($opt_reverse xor $contigmatched) {
        $Nletters_output += $linelength;
        print OUT "$_\n";
    }
}
print STDERR "$Nseqs FASTA sequences in $Nlines lines, total length $Nletters bp\n";
print STDERR "$Nnames unique names, and these matched $Nmatched FASTA sequences\n";
print STDERR "Output $Nletters_output bp\n";
print STDERR "Quit after seeing all $Nnames_seen names\n" if $opt_quitonseen;

