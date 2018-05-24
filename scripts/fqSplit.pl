#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;

my $o_nfiles;
my $o_prefix;
my %filenames;
my %files;
GetOptions(
    "nfiles=i" => \$o_nfiles,
    "prefix=s" => \$o_prefix,
) or die "invalid option: $!";

die "USAGE: $0 --nfiles INT --prefix STRING inputfile1.fastq  (sorry does not currently work for compressed fastq input or output)" if not $o_nfiles or not $o_prefix;

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
