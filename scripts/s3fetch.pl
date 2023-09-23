#!/usr/bin/env perl

# each line is a URL beginning with
# https://s3-us-west-2.amazonaws.com/REMAINDER convert these into
# s3://REMAINDER. So
#
# https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifi/m64076_210310_104300.hifi_reads.fastq.gz
#
# becomes
#
#                               s3://human-pangenomics/T2T/scratch/HG002/sequencing/hifi/m64076_210310_104300.hifi_reads.fastq.gz
#
# Generate and execute aws command lines to download the result.
#
#     aws s3 cp s3://REMAINDER . --no-sign-request
#

use strict;
use warnings;
use feature 'say';

qx/which aws/ or die "could not find aws, module load awscli/1.29.52";

while (<>) {
    chomp;
    (my $s3url = $_) =~ s,https://s3-us-west-2\.amazonaws\.com/,s3://,;
    my $cmd = "aws s3 cp $s3url . --no-sign-request";
    say STDOUT "\n$cmd\n";
    system($cmd);
}
