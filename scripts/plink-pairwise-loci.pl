#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;

sub mysum($) {
    my ($s, $sum) = (@_, 0);
    for (my $i = 0; $i < @$s; ++$i) {
        $sum += $s->[$i];
    }
    return $sum;
}

sub mysum_pairwise_product($$) {
    my ($s, $t, $sum) = (@_, 0);
    die "must be same length" if @$s != @$t;
    for (my $i = 0; $i < @$s; ++$i) {
        $sum += ($s->[$i] * $t->[$i]);
    }
    return $sum;
}

my $o_ped;
my $o_mibs;
my $o_help;

GetOptions("ped:s"  => \$o_ped,
           "mibs:s" => \$o_mibs,
           "help"   => \$o_help) or die "unknown option: $!";

my $usage = "
USAGE:   $0 [ --help ] --ped file.ped --mibs file.mibs

Scan file.ped and write matrix to STDOUT counting the number of complete
pairwise comparisons, with no missing alleles, between each genotype described
in file.ped.  The file.mibs will be used to verify the format of the output but
will not otherwise be used.

NOTE: file.mibs must be derived from file.ped.

A locus with missing alleles is indicated by '0 0' in the PED file.  A message
is produced if only one allele is coded as missing, and the entire locus is
coded as missing.

";
die "$usage" if $o_help or ! $o_ped or ! $o_mibs;

open (my $fd_ped, "<", $o_ped) or die "cannot open PED file '$o_ped': $!";
open (my $fd_mibs, "<", $o_mibs) or die "cannot open MIBS file '$o_mibs': $!";

my %PED;
my %ID;

my $n_loci;
my ($missing, $missing1) = (0, 0);

sub collapse_loci(@) {
    my @a = @_;
    my @l;
    for (my $i = 0; $i < scalar(@a); $i += 2) {
        if ($a[$i] eq '0' and $a[$i + 1] eq '0') {
            ++$missing;
            push @l, 0;
        } elsif ($a[$i] ne '0' and $a[$i + 1] ne '0') {
            push @l, 1;
        } else {
            ++$missing1;
            say STDERR "line $. has 1 of two alleles missing starting at $i";
        }
    }
    $n_loci = scalar(@l) if ! $n_loci;
    die "mismatch at $o_ped line $. between \$n_loci ($n_loci) and number of loci (".scalar(@l).")" if $n_loci != scalar(@l);
    #say STDERR "line $. has ".scalar(@l)." loci of which $m loci are missing both alleles and $m1 are missing one";
    #print Dumper \@l;
    return @l;
}

my $n_ped = 0;
while (<$fd_ped>) {
    s/\R//g;  # used instead of chomp; matches any newline character \n or \r\n
    my ($f, $id, $p1, $p2, $g, $ph, @a) = split /\s+/;
    #print Dumper \$f, \$id, \$p1, \$p2, \$g, \$ph, \@a;
    $PED{$id} = [ collapse_loci(@a) ];
    $ID{$n_ped++} = $id;
}
die "mismatch between \$n_ped ($n_ped) and number of keys in \%ID (".scalar(keys(%ID)).")" if $n_ped != scalar(keys(%ID));

say STDERR "Number of ids read from $o_ped: $n_ped";
say STDERR "Number of loci/id inferred from $o_ped: $n_loci";
say STDERR "Number of loci * number of ids: ".($n_loci * $n_ped);
say STDERR "Number of genotypes missing both alleles: $missing";
say STDERR "Number of genotypes missing one allele: $missing1";

my $n_mibs = 0;
my $n_complete_pairs = 0;
while (<$fd_mibs>) {
    s/\R//g;  # used instead of chomp; matches any newline character \n or \r\n
    my @ibs = split /\s+/;
    die "$o_mibs : $. : empty line" if ! scalar(@ibs);
    my $l_i = $PED{$ID{$n_mibs}};
    my $i_present = mysum($l_i);
    die "MIBS line $. : column $n_mibs is not 1" if $ibs[$n_mibs] != 1;
    die "MIBS line $. : number of columns (".scalar(@ibs).") does not match \$n_ped ($n_ped)" if scalar(@ibs) != $n_ped;
    my @ij;
    for (my $j = 0; $j < scalar(@ibs); ++$j) {
        my $l_j = $PED{$ID{$j}};
        my $j_present = mysum($l_j);
        my $ij_present = mysum_pairwise_product($l_i, $l_j);
        ++$n_complete_pairs if $ij_present == $n_loci;
        die "inconsistent \$i_present ($i_present) or \$j_present ($j_present) or \$ij_present ($ij_present)" if $ij_present > $i_present or $ij_present > $j_present;
        push @ij, $ij_present;
    }
    say STDOUT join(' ', @ij);
    $n_mibs++;
}
say STDERR "Number of rows read from $o_mibs : $n_mibs";
say STDERR "Number of pairwise sets : ".$n_mibs * $n_ped;
say STDERR "Number of complete pairwise sets : $n_complete_pairs";

