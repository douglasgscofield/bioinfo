#!/usr/bin/perl

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use List::Util qw(min max);

sub median { 
	my $ref = shift;
	my @a = sort { $a <=> $b } @{$ref};
	return ($a[$#a/2] + $a[@a/2]) / 2;
}
	
my $infile = $ARGV[0];

open(INFILE, "gzip -f -c -d ${infile} |") or die "couldn't open $infile: $!";
my $fasta_out = Bio::SeqIO->new(-fh => \*STDOUT, -format => "fasta");

my @stats = ();

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
	my $desc = "length:$length,n_reads:$nreads,median_coverage:$medcov";
	my $seqobj = Bio::PrimarySeq->new(-id => $id, -desc => $desc, -seq => $sequence);
	$fasta_out->write_seq($seqobj);
	push @stats, $length;
}

print STDERR "n_seq=" . scalar(@stats) . " min=" . min(@stats) . " max=" . max(@stats) . " median=" . median(\@stats) . "\n";
