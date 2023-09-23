#!/usr/bin/env perl

# report duplicates by digesting .md5 files and comparing MD5 checksums
# requires default modules and Sort::Naturally

use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use Sort::Naturally;

my $USAGE = "
USAGE:   $0  [ --single-line ]  [ --noinclude-md5 ]  [ --noinclude-count ]  [ other options ]  file1.md5  [ file2.md5 ... ]

Scan .md5 lines -- the output of md5sum -- for duplicated MD5 checksums and
report the files having duplicated checksums.  Multiple .md5 files can be
provided on the command line, or on standard input. Output is to standard
output, with informational messages on standard error.

The default output format is multiline, with each set of duplicated files on
separate lines and each block of duplicates separated by one line of ----

    MD5-checksum
    count-of-duplicates
    file1
    file2
    [ ... ]
    ----
    MD5-checksum
    count-of-duplicates
    file1
    file2
    [ ... ]
    ----

OPTIONS

    --single-line      Single-line format: each line contains one set of
                       files sharing the same MD5 checksum, tab separated:

         MD5-checksum <tab> count-of-duplicates <tab> file1 <tab> file2 [ <tab> file3 ... ]

    --noinclude-md5    Do not include the MD5-checksum in the output
    --noinclude-count  Do not include the count of duplicated files in the output

    --noquiet          Produce informational messages
    --help             This help
";

my $o_single_line = 0;
my $o_include_md5 = 1;
my $o_include_count = 1;
my $o_help = 0;
my $o_quiet = 1;

GetOptions("single-line"     => \$o_single_line,
           "include-md5!"    => \$o_include_md5,
           "include-count!"  => \$o_include_count,
           "quiet"           => \$o_quiet,
           "help"            => \$o_help)
or die $USAGE;
if ($o_help) { say $USAGE; exit 0; }

my %M;
my $n_duplicated_md5 = 0;
my $n_duplicated_files = 0;

while (<>) {
    chomp;
    my ($m, $f) = split / +/, $_, 2;
    $f =~ s/^\s+//g;
    $f =~ s/\s+$//g;
    if ($M{$m}) { push @{ $M{$m} }, $f } else { $M{$m} = [ $f ]; }
    say STDERR "processed line $.: '$f'" if ! $o_quiet && $. % 1000 == 0;
}
say STDERR "processed $. lines" if ! $o_quiet;
foreach (sort { ncmp($M{$a}->[0], $M{$b}->[0]) } keys %M) {
    my $n = scalar @{ $M{$_} };
    next if $n == 1;
    ++$n_duplicated_md5;
    $n_duplicated_files += $n;
    if ($o_single_line) {
        print STDOUT "$_\t" if $o_include_md5;
        print STDOUT "$n\t" if $o_include_count;
        print STDOUT "'".join("'\t'", @{ $M{$_} })."'\n";
    } else {
        print STDOUT "$_\n" if $o_include_md5;
        print STDOUT "$n\n" if $o_include_count;
        print STDOUT "'".join("'\n'", @{ $M{$_} })."'\n----\n";
    }
}
if (! $o_quiet) {
    say STDERR "$n_duplicated_md5 duplicated MD5 hashes";
    say STDERR "$n_duplicated_files duplicated files";
    say STDERR sprintf("%.4f files per duplicated MD5 hash", ($n_duplicated_files/$n_duplicated_md5)) if $n_duplicated_md5;
}
