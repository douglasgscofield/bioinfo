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
    print STDERR "
USAGE:  fermiExtractContigs.pl [ options ] fermi.p5.fq.gz

Convert Mummer SNP/indel output as produced by the 'show-snps -T' command to
a pseudo-VCF format.  Indels are collapsed and the first reference base is
inserted at the start of indels as required for VCF, but not all VCF columns
are yet produced.

OPTIONS:

    --minlength    INT     Minimum length of a contig to keep it
    --minreads     INT     Minimum number of reads forming a contig to keep it
    --mincoverage  INT     Minimum median coverage across a contig to keep it
    -h|--help              This help output
    fermi.p5.fq.gz         Fermi-format scaffold output file
";
    exit 1;
}

my $o_help = 0;
my $minlength = 0;
my $n_minlength = 0;
my $minreads = 0;
my $n_minreads = 0;
my $mincoverage = 0;
my $n_mincoverage = 0;

GetOptions('minlength=i' => \$minlength,
           'minreads=i' => \$minreads,
           'mincoverage=i' => \$mincoverage,
           'help' => \$o_help) or usage("unrecognized option");
usage("") if $o_help;

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
	my ($id, $length, $nreads) = split /\t/;
	$id = substr($id, 1); # drop '@'
	my $sequence = <INFILE>;
	chomp $sequence;
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
	my $desc = "length:$length,n_reads:$nreads,median_coverage:$medcov";
	my $seqobj = Bio::PrimarySeq->new(-id => $id, -desc => $desc, -seq => $sequence);
	$fasta_out->write_seq($seqobj);
	push @stats, $length;
}

print STDERR "KEPT: n_seq=", scalar(@stats), " total length=", sum(@stats), 
             " min=", min(@stats), " max=", max(@stats), " median=", median(\@stats), "\n",
             "REJECTED: n_seq=", scalar(@rejects),  " total length=", sum(@rejects),
             " min=", min(@rejects), " max=", max(@rejects), 
             " median=", median(\@rejects), "\n",
             "N sequences rejected for length < $minlength: ", $n_minlength, "\n",
             "N sequences rejected for N reads < $minreads: ", $n_minreads, "\n",
             "N sequences rejected for median coverage < $mincoverage: ", $n_mincoverage, "\n";
