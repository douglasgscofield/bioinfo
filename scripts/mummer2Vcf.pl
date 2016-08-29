#!/usr/bin/env perl

# Changelog
#
# 2016-04-19: Liangjiao Xue <lxue@uga.edu> caught a bug where the alignment
#             direction was not honoured when adding indel sequence.  Thanks!

use warnings;
use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

sub usage($) {
    my $msg = shift;
    print STDERR "ERROR: $msg\n" if $msg;
    print STDERR "
USAGE:  mummer2Vcf.pl [ options ] mummer.snps

Convert Mummer SNP/indel output as produced by the 'show-snps -T' command to
a pseudo-VCF format.  Indels are collapsed and the first reference base is
inserted at the start of indels as required for VCF, but not all VCF columns
are yet produced.

If indels are to be output, the reference sequence(s) against which Mummer
called SNPs and indels must be either available as the file named at the start
of the first line in the mummer.snps file, or via the -f/--fasta option.

OPTIONS:

    -f|--fasta  FILE.fa     Fasta file for reference, in same order as
                            the reference in the mummer.snps file
    -t|--type   SNP|INDEL   Restrict the output to just SNPs or indels;
                            if the output is just SNPs, no reference
                            file is required
    -h|--help               This help output
    mummer.snps             SNPs/indels as output by the Mummer command
                            'show-snps -T'
";
    exit 1;
}

my $indel = "";
my $indel_type = "";
my $indel_ref = "";
my $indel_ref_start = 0;
my $indel_query = "";
my $indel_query_start = 0;
my $fasta_file = "";
my $o_type = "";
my $o_snpEffect = 0;
my $o_help = 0;
my $ref_fasta_io;
my $ref_fasta_seq;
my $align_direction = 1;

GetOptions('fasta=s' => \$fasta_file,
           'type=s' => \$o_type,
           'snpEffect' => \$o_snpEffect,
           'help' => \$o_help) or usage("");
usage("") if $o_help;
usage("--type must be SNP or INDEL") if $o_type and $o_type ne "INDEL" and $o_type ne "SNP";

if ($o_type ne "SNP") {
    # only open the Fasta file if we're including indels in the output
    if ($fasta_file eq "") {
        my $firstline = <>;
        ($fasta_file, undef) = split /\s+/, $firstline;
    }
    usage("must specify Fasta reference with -f/--fasta <file>") if $fasta_file eq "";
    $ref_fasta_io = Bio::SeqIO->new(-file => "<$fasta_file", -format => "fasta");
}

sub find_ref_fasta_seq($) {
    return if $o_type eq "SNP";
    # read fasta sequences until we load the requested one
    my $seq_name_to_find = shift;
    # if we return below, then we are already at the sequence of interest
    if (defined($ref_fasta_seq) and $ref_fasta_seq->id eq $seq_name_to_find) {
        #print STDERR "find_ref_fasta_seq: we are already at the reference sequence $seq_name_to_find\n";
        return;
    }
    while (my $s = $ref_fasta_io->next_seq()) {
        if ($s->id eq $seq_name_to_find) {
            $ref_fasta_seq = $s;
            #print STDERR "find_ref_fasta_seq: found new reference sequence $seq_name_to_find\n";
            return;
        }
    }
    die("couldn't find reference $seq_name_to_find, are references for SNPs/indels and sequences in --fasta file in the same order?");
}
sub complete_indel() {
    if ($o_type ne "SNP") {
        # fetch the ref base from the reference sequence
        find_ref_fasta_seq($indel_ref);
        my $ref_base = $ref_fasta_seq->subseq($indel_ref_start - 1, $indel_ref_start - 1);
        $indel = $ref_base . $indel;
        my $alt = $ref_base;
        if ($indel_type eq "INSERTION") {
            #print STDOUT join("\t", $indel_ref, $indel_ref_start, $alt, $indel, $indel_type), "\n";
            print STDOUT join("\t", $indel_ref, $indel_ref_start, $alt, $indel, "INDEL"), "\n";
        } elsif ($indel_type eq "DELETION") {
            #print STDOUT join("\t", $indel_ref, $indel_ref_start, $indel, $alt, $indel_type), "\n";
            print STDOUT join("\t", $indel_ref, $indel_ref_start, $indel, $alt, "INDEL"), "\n";
        } else {
            die("complete_indel: oops");
        }
    }
    $indel = "";
    $indel_type = "";
    $indel_ref = "";
    $indel_ref_start = 0;
    $indel_query = "";
    $indel_query_start = 0;
}
sub complete_snp($$$$) {
    if ($o_type ne "INDEL") {
        my ($ref, $pos, $orig, $alt) = @_;
        print STDOUT join("\t", $ref, $pos, $orig, $alt, "SNP"), "\n";
    }
}

while (<>) {
    last if (/^NUCMER/);  # skip all lines up to NUCMER, as produced in proper Mummer output
}
if ($o_snpEffect) {
    # print snpEffect header, first 4 (strictly middle 3) are required
    print STDOUT join("\t", "Reference", "Position", "ReferenceAllele", "SNPAllele", "Type"), "\n";
} else {
    # print VCF header
}
# Note: all --type processing is handled in the complete_{indel,snp}() subroutines.  We
# still track indels even when we don't output them.
while (<>) {
    next if /(^NUCMER|^\s*$|^\[)/;  # skip headers and blank lines
    my @line = split;
    if ($indel and (($line[1] eq "." and $line[0] != $indel_ref_start) 
                    or ($line[2] eq "." and $line[3] != $indel_query_start))) {
        # we have another indel but a new one
        complete_indel();
    }
    if( $line[1] eq "." || $line[2] eq ".") {
        if ($indel eq "") {
            # a new indel, track it
            $indel_ref = $line[10];
            $indel_ref_start = $line[0];
            $indel_query = $line[11];
            $indel_query_start = $line[3];
            $indel_type = $line[1] eq "." ? "INSERTION" : "DELETION";
            $align_direction = $line[9];
        }
        if ($line[1] eq ".") { # insertion in query
            if($align_direction == 1){
            	 $indel .= $line[2];
            } else{
            	 $indel = $line[2].$indel;
            }
        } elsif ($line[2] eq ".") { # deletion in query
            $indel .= $line[1];
        } else {
            die("oops");
        }
    } else {
        # a snp
        complete_indel() if $indel;
        complete_snp($line[10], $line[0], $line[1], $line[2]);
    }
}
complete_indel() if $indel;

