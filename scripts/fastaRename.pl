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

my $o_infile = "";
my $o_outfile = "";
my $o_mapfile = "";
my $o_prefix = "";
my $o_nozeropad = 0;
my $o_description = "shift";
my @description_actions = qw/shift keep strip/;
my $o_width = 0;
my $o_stdio;
my $o_help;

my $usage = "
NAME

    $script - Rename fasta sequences using simple prefix-number scheme

SYNOPSIS

    Uses only standard Perl modules.

    Fasta sequences are renamed to PREFIX###### where ###### is a zero-padded
    sequence number, starting with 1.  PREFIX is specified with the --prefix
    option.  The number of digits is specified with --width or, if it is not
    specified and the input is a file and not stdin, the number of sequences
    in the input file is counted and the value is calculated as the minimum
    number of zero-padded digits required to number all sequences.

    If --no-zero-pad is specified, then --width is not required, the sequences
    are not pre-counted, and the sequence number is suffixed to PREFIX without
    padding.

    The --prefix option must be specified.

        $script --prefix Trebouxia --in fasta.fa --out fasta.renamed.fa --map fasta.renamed.map.txt

        $script --prefix metascaffold -w 10 --in fasta.fa --out fasta.renamed.fa

        $script --prefix seq --no-zero-pad - < fasta.fa > fasta.renamed.fa

        cat fasta.fa | $script --prefix Trebouxia --width 5 - --out renamed.fa

OPTIONS

    -p PREFIX, --prefix PREFIX   REQUIRED. Use PREFIX in the names of the renamed
                                 sequences.

    -w WIDTH, --width WIDTH      Use WIDTH zero-padded digits to number the
                                 renamed sequences. If WIDTH is not specified and
                                 the input is a file, then the sequences are 
                                 pre-counted from the file and WIDTH is calculated
                                 from the sequence count; this requires a bit of
                                 time to complete.

    -n, --no-zero-pad            Do not zero-pad sequence numbers.  Sequence numbers
                                 are not pre-counted, and --width is not allowed.

    -d ACTION, --description ACTION
                                 ACTION describes how to handle the current sequence
                                 '>name description'. A fasta description is
                                 everything on the name line after the first
                                 whitespace. The options are case-insensitive.

                          Default -> SHIFT : Shift 'name description' to become the
                                             description for the new sequence. This
                                             keeps the old name next to the new name:
                                             '>newname###### name description'
                                     KEEP  : Keep the description as the description
                                             for the renamed sequence. Discard the
                                             old name:
                                             '>newname###### description'
                                     STRIP : Strip all description information from
                                             the renamed sequence:
                                             '>newname######'
                                 
    -m FILE, --map FILE          Write a sequence name map to FILE containing lines
                                 of the format:
                                 'newname###### TAB name'

    -                            Read FASTA sequences from stdin and/or write
                                 extracted sequences to stdout
    -i FILE, --in FILE           File for input FASTA sequences
    -o FILE, --out FILE          File to output FASTA sequences

    -?, --help                   This help message
";

sub print_usage_and_exit($) {
    my $code = shift;
    say $usage;
    exit $code;
}

print_usage_and_exit(1) if ! scalar @ARGV;

GetOptions(
    "" => \$o_stdio, 
    "i|in=s" => \$o_infile,
    "o|out=s" => \$o_outfile,
    "m|map=s" => \$o_mapfile,
    "p|prefix=s" => \$o_prefix,
    "n|nozeropad|no-zero-pad" => \$o_nozeropad,
    "d|description" => \$o_description,
    "w|width=i" => \$o_width,
    "help|?" => \$o_help,
) or print_usage_and_exit(1);

print_usage_and_exit(0) if $o_help;

$o_infile  = "/dev/stdin"  if ! $o_infile and $o_stdio;
$o_outfile = "/dev/stdout" if ! $o_outfile and $o_stdio;

print_usage_and_exit(1) if ! $o_infile or ! $o_outfile;
print_usage_and_exit(1) if $o_width and $o_nozeropad;

if ($o_infile eq "/dev/stdin" and ! $o_nozeropad and $o_width == 0) {  # quickly count the number of sequences
    say STDERR "*** must specify --width if reading from stdin without --no-zero-pad";
    print_usage_and_exit(1);
}

@description_actions = grep { $o_description =~ m/^$_$/i } @description_actions;
die "unknown --description action" if scalar @description_actions != 1;
$o_description = lc($description_actions[0]);

if ($o_width == 0 and ! $o_nozeropad) {
    say STDERR "$script: counting sequences in $o_infile ...";
    my $n = qx/grep -c '^>' "$o_infile"/;
    chomp $n;
    die "Could not count fasta sequences in '$o_infile': the file is not found, 'grep' is not found, or something else: $!" if !$n;
    $o_width = int(log($n)/log(10)) + 1;
    say STDERR "$script: found $n sequences, calculating $o_width digits for WIDTH";
}
open(my $infd,  "<$o_infile")  or die("couldn't open fasta input: $!");
open(my $outfd, ">$o_outfile") or die("couldn't open fasta output: $!");

my $mapfd;
if ($o_mapfile) {
    open($mapfd, ">$o_mapfile") or die("couldn't create --map file '$o_mapfile': $!");
}

my $newname_template = $o_nozeropad ?  "${o_prefix}%d" : "${o_prefix}%0${o_width}d";
my $n_sequence = 0;

while (<$infd>) {
    chomp;
    if (/^>/) {
        ++$n_sequence;
        my $header = substr($_, 1);
        my ($name, $description) = split /\s+/, $header, 2;
        my $newname = sprintf($newname_template, $n_sequence);
        my $newheader;
        if ($o_description eq 'shift') {
            $newheader = "${newname} ${header}";
        } elsif ($o_description eq 'keep' and $description) {
            $newheader = "${newname} ${description}";
        } else {
            $newheader = "${newname}";
        }
        say {$outfd} ">$newheader";
        say {$mapfd} "$newname\t$name" if $o_mapfile;
    } else {
        say {$outfd} "$_";
    }
}
$infd->close();
$outfd->close();
$mapfd->close() if $o_mapfile;
say STDERR "$script: sequences renamed : $n_sequence";
say STDERR "$script: rename template   : '$newname_template'";
say STDERR "$script: input             : $o_infile";
say STDERR "$script: output            : $o_outfile";
say STDERR "$script: names map         : $o_mapfile" if $o_mapfile;

