#!/usr/bin/env perl

# Convert 'exonerate --model protein2genome --showfeaturegff true' output to CDS DNA sequence and translated protein.
#
# arg 1: Genome fasta containing the genome DNA sequences referred to by the GFF
# arg 2: GFF produced by 'exonerate --model protein2genome --showfeaturegff true'
# arg 3: prefix for the two output files:
#        1. concatenated DNA sequences from segments marked 'cds' in the GFF, named prefix.cds.fasta
#        2. translated version of this sequence, named prefix.pep.fasta
#
# Requires BioPerl.  BioPerl is used for indexing the genome fasta, extracting
# CDS segments from the genome fasta, doing reverse-complementing if required,
# translating the sequence to proteins, and writing the output files.
#
# The gff output by exonerate is *not* gff like output by maker.  There's no
# mRNA line and there are extra lines marking introns, splice sites, etc.  This
# script is specific to exonerate's output (as of version 2.4.0)
#
# Derived from a script posted by Federico Giorgi at https://biostars.org/p/46281 with some modifications
#   - moved code into start_gene() and wrapup_gene()
#   - each sequence is tagged with gene_<gene # from exonerate output>_query_<query name>
#   - fixed last-gene-not-output bug


use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

if (! $ARGV[0]) {
    print STDERR "$0 file.fasta file.gff3 output-prefix\n";
    exit 1;
}

my $file_fasta = $ARGV[0];
my $file_gff = $ARGV[1];
my $output_prefix = $ARGV[2];

$| = 1;    # Flush output
my $outfile_cds = Bio::SeqIO->new( -format => 'fasta', -file => ">${output_prefix}.cds.fasta" );
my $outfile_pep = Bio::SeqIO->new( -format => 'fasta', -file => ">${output_prefix}.pep.fasta" );

### First, index the genome
my $db = Bio::DB::Fasta->new($file_fasta);
print STDERR "Genome fasta parsed\n";

open GFF, "<$file_gff" or die $!;

### Second, parse the GFF3
my @mRNA;
my $mRNA_name;
my $frame;

sub start_gene(@) {
    my @array = @_;
    # Now initialize the next mRNA
    my @attrs = split( " ; ", $array[8] );
    $attrs[0] =~ s/gene_id //;
    $attrs[1] =~ s/sequence //;
    $attrs[2] =~ s/gene_orientation //;
    $mRNA_name = "gene_$attrs[0]_query_$attrs[1]";
    $frame = $attrs[2];
    @mRNA = (); # Empty the mRNA
}

sub wrapup_gene {
    my @array = @_;
    # Collect CDSs and extract sequence of the previous mRNA
    my $mRNA_seq;
    foreach my $coord (@mRNA) {
        my @cds_coord = split( " ", $coord );
        my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
        $mRNA_seq .= $cds_seq;
    }

    my $output_nucleotide = Bio::Seq->new(
        -seq        => $mRNA_seq,
        -id         => $mRNA_name,
        -display_id => $mRNA_name,
        -alphabet   => 'dna',
    );
    if ($frame eq '-') {
        $output_nucleotide = $output_nucleotide->revcom();
    }
    my $output_protein = $output_nucleotide->translate();
    $outfile_cds->write_seq($output_nucleotide);
    $outfile_pep->write_seq($output_protein);

    # Now initialize the next mRNA
    start_gene(@array) if @array;
    #my @attrs = split( ";", $array[8] );
    #$attrs[0] =~ s/gene_id //;
    #$attrs[1] =~ s/sequence //;
    #$attrs[2] =~ s/gene_orientation //;
    #$mRNA_name = $attrs[1];
    #$frame=$array[2];
    #@mRNA = (); # Empty the mRNA
}

while ( my $line = <GFF> ) {
    next if $line =~ /^#/ or $line =~ /completed exonerate analysis/;
    chomp $line;
    my @array = split( "\t", $line );
    my $type = $array[2];
    next if ( $type eq 'splice5' or $type eq 'splice3' or $type eq 'intron' or $type eq 'UTR' or $type eq 'similarity' );

    if ( ( $type eq 'gene' and @mRNA ) ) {
        wrapup_gene(@array);
    }
    elsif ( $type eq 'gene' ) {    # First mRNA
        start_gene(@array);
    }
    elsif ( $type eq 'cds' ) {
        my $cds_coord = $array[0] . " " . $array[3] . " " . $array[4];
        push( @mRNA, $cds_coord );
    }
}
wrapup_gene();

close GFF;
