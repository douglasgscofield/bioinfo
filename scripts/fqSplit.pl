#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;

my $o_nfiles;
my $o_prefix;
my $o_help;
my %filenames;
my %files;
GetOptions(
    "nfiles=i" => \$o_nfiles,
    "prefix=s" => \$o_prefix,
    "help"     => \$o_help,
) or die "invalid option: $!";

die "
    USAGE :  $0 [ --help ] --nfiles INT --prefix STRING reads.fastq

With compressed input/output:

    USAGE : zcat reads.fastq.gz | $0 --nfiles INT --prefix STRING && gzip -v STRING*.fastq


Read reads.fastq and split it into INT different files, each filename starting
with STRING, followed by the padded file number (1 to INT) and then ending with
the suffix '.fastq'.

Sorry, the script does not currently support compressed input or output, so
some variation of the second form with the pipe and ending with '&& gzip ...'
can be used to uncompress reads for the script and then compress the output.

" if $o_help or (not $o_nfiles and not $o_prefix);

die "Must provide values for both --nfiles INT and --prefix STRING" if not $o_nfiles or not $o_prefix;

my $padding = int(log($o_nfiles) / log(10)) + 1;

for my $i (1 .. $o_nfiles) {
    $filenames{$i} = sprintf("${o_prefix}%.${padding}d.fastq", $i);
    say STDERR "filename $i is $filenames{$i}";
    open($files{$i}, ">", $filenames{$i}) or die "could not create file $filenames{$i}: $!";
}

my $n_reads = 0;
while (<>) {
    my $read_file = ( ($. - 1) / 4  % $o_nfiles) + 1;
    #say "read the $.-th line of input, writing to the $read_file-th file";
    my $l1 = $_;
    my $l2 = <>;
    my $l3 = <>;
    my $l4 = <>;
    print { $files{$read_file} } $l1;
    print { $files{$read_file} } $l2;
    print { $files{$read_file} } $l3;
    print { $files{$read_file} } $l4;
    ++$n_reads;
}
$files{$_}->close() foreach keys %files;

say STDERR "split $n_reads reads into $o_nfiles files";
