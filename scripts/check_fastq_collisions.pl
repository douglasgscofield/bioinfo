#!/usr/bin/env perl

# Receive a list of filenames on stdin, check for collisions among *.fastq.{gz,bz2}
# filenames and among dir/*.fastq.{gz,bz2} filenames.  If there is a collision, or
# a duplicate path given on input, 'exit 1' else 'exit 0'
#
# With -v/--verbose, print the full paths of the duplicates
#
# It would be nice if collisions excluding the .gz/.bz2 were detected

use strict;
use warnings;
use File::Spec;
use Getopt::Long;

my $o_verbose = 0;
my $dups_seen = 0;

my %fullpaths;
my %files;
my %dirfiles;
my %map;

GetOptions("verbose" => \$o_verbose) or die;

while (my $path = <>) {
    chomp $path;
    ++$fullpaths{$path};
    next if $path !~ /\.fastq\.(gz|bz2)$/;  # not *.fastq.{gz,bz2}
    my (undef, $directories, $file) = File::Spec->splitpath($path);
    $directories =~ s/\/$//;  # remove trailing / if present
    ++$files{$file};
    push @{$map{$file}}, $path;
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
check_duplicates(\%fullpaths, "full input path");
check_duplicates(\%files,     "filename");
check_duplicates(\%dirfiles,  "directory/filename");

exit 1 if $dups_seen;
