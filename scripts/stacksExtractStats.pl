#!/usr/bin/env perl

# Extract sample, stack and stack coverage information from the Stacks log file.
# For the following input:
# 
#     Identifying unique stacks; file   1 of  40 [XYZ_H0122_GATTAC.1.cat]
#     /xxxx/yyyy/private/GBS_working_data/stacks-1.20/ustacks -t fastq -f ./XYZ/concat/XYZ_H0122_GATTAC.1.cat.fq -o ./XYZ_read1 -i 1 -m 3 -M 3 -p 15 -d -r 2>&1
#     Min depth of coverage to create a stack: 3
#     Max distance allowed between stacks: 3
#     Max distance allowed to align secondary reads: 5
#     Max number of stacks allowed per de novo locus: 3
#     Deleveraging algorithm: enabled
#     Removal algorithm: enabled
#     Model type: SNP
#     Alpha significance level for model: 0.05
#     Parsing ./XYZ/concat/XYZ_H0122_GATTAC.1.cat.fq
#     Loaded 2800 RAD-Tags; inserted 2423 elements into the RAD-Tags hash map.
#       0 reads contained uncalled nucleotides that were modified.
#       Mean coverage depth is 9; Std Dev: 18.0875 Max: 87
#     Coverage mean: 9; stdev: 18.0875
#     Deleveraging trigger: 27; Removal trigger: 45
#     Calculating distance for removing repetitive stacks.
#       Distance allowed between stacks: 1
#       Using a k-mer length of 71
#       Number of kmers per sequence: 74
#       Miniumum number of k-mers to define a match: 3
#     Removing repetitive stacks.
#       Removed 3 stacks.
#       38 stacks remain for merging.
#     Calculating distance between stacks...
#       Distance allowed between stacks: 3
#       Using a k-mer length of 35
#       Number of kmers per sequence: 110
#       Miniumum number of k-mers to define a match: 5
#     Merging stacks, maximum allowed distance: 3 nucleotide(s)
#       38 stacks merged into 36 stacks; deleveraged 0 stacks; removed 0 stacks.
#       Mean merged coverage depth is 9.41667; Std Dev: 20.8037; Max: 96
#     Merging remainder radtags
#       2461 remainder sequences left to merge.
#       Distance allowed between stacks: 5
#       Using a k-mer length of 23
#       Number of kmers per sequence: 122
#       Miniumum number of k-mers to define a match: 7
#       Matched 29 remainder reads; unable to match 2432 remainder reads.
#     Number of utilized reads: 368
#     Writing loci, SNPs, and alleles to './XYZ_read1/'...
#       Refetching sequencing IDs from ./XYZ/concat/XYZ_H0122_GATTAC.1.cat.fq... read 2800 sequence IDs.
#     done.
#    
# The command
#
#     ./stacksExtractStats.pl ex.log
#
# produces
#
#     samplename	nsample	totsamples	nstacks	coverage
#     XYZ_H0122_GATTAC.1	1	40	36	9.41667
# 

use strict;
use warnings;

print "samplename\tnsample\ttotsamples\tnstacks\tcoverage\n";

my $totsamples;

while (<>) {
    my ($nsample, $samplename, $nstacks, $coverage);
    s/\R//g;  # general line ending chomper: both \n and \r\n
    # the regexp matches the first sample-specific line in the log for each sample
    if (/^Identifying unique stacks; file +([0-9]+) of +([0-9]+) \[(.*)\.cat\]$/) {
        $nsample = $1;
        $samplename = $3;
        if (! $totsamples) {
            $totsamples = $2;
        } elsif ($2 != $totsamples) {
            die "sample total at sample number $nsample, name $samplename: $2 vs. $totsamples";
        }
        while (<>) {
            s/\R//g;  # general line ending chomper: both \n and \r\n
            if (/^Error: Unable to form any stacks/) {
                $nstacks = 0;
                $coverage = 0;
            } elsif (/^  [0-9]+ stacks merged into ([0-9]+) stacks;/) {
                $nstacks = $1;
                $_ = <>;
                s/\R//g;  # general line ending chomper: both \n and \r\n
                s/^  Mean merged coverage depth is ([0-9.]+);.*$/$1/g;
                $coverage = $_;
            } else {  # some other line, skip and read next
                next;
            }
            print "$samplename\t$nsample\t$totsamples\t$nstacks\t$coverage\n";
            last;
        }
    }
}
