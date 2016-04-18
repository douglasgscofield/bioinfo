#!/usr/bin/env perl

# Receive a list of filenames on stdin, check for collisions among
# *.fastq.{gz,bz2} filenames and among dir/*.fastq.{gz,bz2} filenames.  If
# there is a collision, or a duplicate path given on input, 'exit 1' else 'exit
# 0'
#
# With -v/--verbose, print the full paths of the duplicates
#
# With --trimsuffix INT, duplicates of filename after INT dot-separated
# suffixes are removed are also checked.  This can be used to check for FastQ
# files that collide except for compression method, such as file.fastq.gz and
# file.fastq.bz2 with --trimsuffix 1, and files that collide except for FastQ
# suffix, such as file.fastq.gz and file.fq.gz with --trimsuffix 2.

use strict;
use warnings;
use File::Spec;
use Getopt::Long;

my $o_verbose = 0;
my $o_trimsuffix = 0;
my $dups_seen = 0;

my %fullpaths;
my %files;
my %files_trim;
my %dirfiles;
my %map;

GetOptions("verbose"      => \$o_verbose,
           "trimsuffix=i" => \$o_trimsuffix) or die "USAGE: cat list-of-filenames | $0 [ -v ] [ --trimsuffix INT ]";

while (my $path = <>) {
    chomp $path;
    ++$fullpaths{$path};
    next if $path !~ /\.(fastq|fq)\.(gz|bz2|xz)$/;  # not *.fastq.{gz,bz2,xz}
    my (undef, $directories, $file) = File::Spec->splitpath($path);
    $directories =~ s/\/$//;  # remove trailing / if present
    ++$files{$file};
    push @{$map{$file}}, $path;
    if ($o_trimsuffix) {
        # split filename by dot, join after removing the last $o_trimsuffix
        my @t = split(/\./, $file);
        if (scalar @t > $o_trimsuffix) {
            pop @t for 1 .. $o_trimsuffix;
            my $filetrim = join('.', @t);
            ++$files_trim{$filetrim};
            push @{$map{$filetrim}}, $path;
        }
    }
    my @dirs = File::Spec->splitdir($directories);
    if (scalar @dirs) {
        my $dirfile = File::Spec->catfile($dirs[$#dirs], $file);
        ++$dirfiles{$dirfile};
        push @{$map{$dirfile}}, $path;
    }
}
sub check_duplicates($$) {
    my ($hash, $tag) = @_;
    for my $k (sort keys %{$hash}) {
        if ($hash->{$k} > 1) {
            print "duplicate $tag x$hash->{$k}: $k\n";
            if ($o_verbose) {
                 print "\t$_\n" foreach (@{$map{$k}});
            }
            ++$dups_seen;
        }
    }
}
check_duplicates(\%fullpaths,  "full input path");
check_duplicates(\%files,      "filename");
check_duplicates(\%files_trim, "filename removing $o_trimsuffix suffixes") if $o_trimsuffix;
check_duplicates(\%dirfiles,   "directory/filename");

exit 1 if $dups_seen;
