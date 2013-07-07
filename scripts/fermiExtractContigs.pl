#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use List::Util qw(min max sum);

sub usage($) {
    my $msg = shift;
    print STDERR "ERROR: $msg\n" if $msg;
    print STDERR <<'__USAGE__';

Usage:  fermiExtractContigs.pl [ options ] file.{fq,mag}.gz

Extract Fasta-format contigs from Fermi's FastQ-like output files.  Fermi is an
overlap assembler developed by Heng Li (https://github.com/lh3/fermi).  

Fermi produces two formats of FastQ-like output: scaftig sequences in files
named *.fq.gz, and unitig sequences in files named *.mag.gz which use an
overlap graph file format called MAG.  This script can read both formats, with
scaftigs being the default input and MAG format interpreted with the --mag or
--mag-verbose options.

Output of this script is annotated, line-wrapped Fasta format, to stdout.  For
scaftig files, the FastQ-like read name line includes two additional
tab-separated fields encoding the length of the scaftig and the number of
non-redundant reads that built the scaftig.  This script includes these in the
Fasta description, and also calculates the median non-redundant read coverage
along the contig by interpreting the encoded coverage in fourth line.  So this:

    @26417937:25351227_0	191	83
    TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTATATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTAAACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCCTGCATGTAGAA
    +
    \"#$%%&'((()**+,-../0123456789:;<=>?@ABCDEFGHIJKLLMNOPQRSTUVWXYZ[\]^_`abcdeffghijklmnoopqrstttttttttttsrqpponmmmlkkjihggfedcba`_^]\[ZYXWVUTSRQPONMLKJIIHGFEDCBA@?>=<;:9876543210//.-,+*)('&&%$#\"

becomes this:

    >26417937:25351227_0 length:191;n_reads:83;median_coverage:44
    TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTAT
    ATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTA
    AACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCC
    TGCATGTAGAA

MAG-format input presents different fields in the read name line.  In MAG, the
second field is the number of non-redundant reads, and the following two fields
list the left and right neighbors of the unitig in the overlap graph structure.
Unitigs may have multiple left or right neighbours, or none, in which case the
neighbour information is a single '.'.  Because unitig length is not encoded in
the header, it is calculated from each unitig.  In the example below, 13 reads
were used to build the unitig, and it has two neighbour unitigs to the left,
one with 77 and one with 76 bases of exact overlap, and four neighbor unitigs
to the right, with 67, 60, 68 and 69 bases of exact overlap.

    @9012595342:2229517731	13	355005708,77;958384802,76;	2817374678,67;3657728097,60;7691722363,68;9232081066,69;
    TGTTTTATTAATAAAATATCTCTCTAACTTGTTTATTTCTGACCTGTTTCAGGTGAATTCGAAATTGCATGGGATTTGATGGATAGTTTAGAGGATATTTGGATGATTTATTTTCATTTTAGTTTCCTAGTTTAGCTGATCTTGGGAATATCTTCCCAGATTATGTAAACAGTTTGATAACTTCCACAGGGAGGTTTTACCCTGTGGAAAACTTATAAATACTTATTAT
    +
    """""""""##########$$$$$%%&&&&&&&&&&&&&&&&&&&&&&&&&'''''''''(((((((((()))))))*******+++++++++++++++++*********)))))))***)))))))())))))))))))))))))))))))(((((((((''''''''''&&&&&&&%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$########"""

Output with --mag is

    >9012595342:2229517731 length:229;n_reads:13;median_coverage:6
    TGTTTTATTAATAAAATATCTCTCTAACTTGTTTATTTCTGACCTGTTTCAGGTGAATTC
    GAAATTGCATGGGATTTGATGGATAGTTTAGAGGATATTTGGATGATTTATTTTCATTTT
    AGTTTCCTAGTTTAGCTGATCTTGGGAATATCTTCCCAGATTATGTAAACAGTTTGATAA
    CTTCCACAGGGAGGTTTTACCCTGTGGAAAACTTATAAATACTTATTAT

Output with --mag-verbose includes the lists of left and right neighbours:

    >9012595342:2229517731 length:229;n_reads:13;median_coverage:6 355005708,77;958384802,76; 2817374678,67;3657728097,60;7691722363,68;9232081066,69;
    TGTTTTATTAATAAAATATCTCTCTAACTTGTTTATTTCTGACCTGTTTCAGGTGAATTC
    GAAATTGCATGGGATTTGATGGATAGTTTAGAGGATATTTGGATGATTTATTTTCATTTT
    AGTTTCCTAGTTTAGCTGATCTTGGGAATATCTTCCCAGATTATGTAAACAGTTTGATAA
    CTTCCACAGGGAGGTTTTACCCTGTGGAAAACTTATAAATACTTATTAT


OPTIONS:

    --mag                 Input is fermi MAG format for unitigs instead of 
                          scaffold FastQ
    --mag-verbose         Input is fermi MAG format, include lists of left and
                          right unitig neighbours in output
    --minlength    INT    Minimum length of a contig to keep it
    --minreads     INT    Minimum number of reads forming a contig to keep it
    --mincoverage  INT    Minimum median coverage across a contig to keep it
    -h|--help             This help output

    file.fq.gz            Fermi-format scaffold file as input
    file.mag.gz           Fermi-format MAG unitig file as input

This script required BioPerl 1.6.1 to compose and write the Fasta-format
output.

__USAGE__

    exit 1;
}

my $o_help = 0;
my $minlength = 0;
my $n_minlength = 0;
my $minreads = 0;
my $n_minreads = 0;
my $mincoverage = 0;
my $n_mincoverage = 0;
my $o_mag = 0;
my $o_mag_verbose = 0;

GetOptions('mag' => \$o_mag,
           'mag-verbose' => \$o_mag_verbose,
           'minlength=i' => \$minlength,
           'minreads=i' => \$minreads,
           'mincoverage=i' => \$mincoverage,
           'help' => \$o_help) or usage("unrecognized option");
usage("") if $o_help;
$o_mag = 1 if $o_mag_verbose;

sub median { 
    my $ref = shift;
    my @a = sort { $a <=> $b } @{$ref};
    if (scalar(@a) % 2 == 0) {
        return (($a[$#a/2] + $a[($#a/2)+1]) / 2);
    } else {
        return ($a[$#a/2]);
    }
}
    
my $infile = $ARGV[0];
if (! $infile) {
    usage("");
}

open(INFILE, "gzip -f -c -d ${infile} |") or die "couldn't open $infile: $!";
my $fasta_out = Bio::SeqIO->new(-fh => \*STDOUT, -format => "fasta");

my @stats = ();
my @rejects = ();


while (<INFILE>) {
    chomp;
    my ($id, $length, $nreads, $left, $right);
    if ($o_mag) {
        ($id, $nreads, $left, $right) = split /\t/;
    } else {
        ($id, $length, $nreads) = split /\t/;
    }
    $id = substr($id, 1); # drop '@'
    my $sequence = <INFILE>;
    chomp $sequence;
    $length = length($sequence) if $o_mag;  # not supplied in header in MAG
    scalar(<INFILE>); # throw away the '+' line
    my $coverage = <INFILE>;
    chomp $coverage;
    my @cov = map { $_ - 33 } unpack('U*', $coverage);
    my $medcov = median(\@cov);
    if ($minlength and $length < $minlength) {
        ++$n_minlength; 
        push @rejects, $length;
        next;
    }
    if ($minreads and $nreads < $minreads) {
        ++$n_minreads;
        push @rejects, $length;
        next;
    }
    if ($mincoverage and $medcov < $mincoverage) {
        ++$n_mincoverage;
        push @rejects, $length;
        next;
    }
    my $desc = "length:$length;n_reads:$nreads;median_coverage:$medcov";
    $desc .= " $left $right" if $o_mag_verbose;
    my $seqobj = Bio::PrimarySeq->new(-id => $id, -desc => $desc, -seq => $sequence);
    $fasta_out->write_seq($seqobj);
    push @stats, $length;
}

print STDERR "KEPT: n_seq=", scalar(@stats), " total length=", sum(@stats), 
             " min=", min(@stats), " max=", max(@stats), " median=", median(\@stats), "\n";
if (scalar(@rejects) > 0) {
    print STDERR "REJECTED: n_seq=", scalar(@rejects),  " total length=", sum(@rejects),
                " min=", min(@rejects), " max=", max(@rejects), 
                " median=", median(\@rejects), "\n",
                "N sequences rejected for length < $minlength: ", $n_minlength, "\n",
                "N sequences rejected for N reads < $minreads: ", $n_minreads, "\n",
                "N sequences rejected for median coverage < $mincoverage: ", $n_mincoverage, "\n";
} else {
    print STDERR "REJECTED: none\n";
}
