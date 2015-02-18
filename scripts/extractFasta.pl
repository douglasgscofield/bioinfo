#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Seq;  # the BioPerl sequence manipulation package
use Bio::SeqIO;  # the BioPerl sequence file manipulation package
use Getopt::Long;

my $o_db;
my $o_entry;
my $o_range;

GetOptions("db=s" => \$o_db,
           "entry=s" => \$o_entry,
           "range=s" => \$o_range) or die "bad option";
die "must supply --db and --entry, optionally --range" if not $o_db or not $o_entry;
print STDERR "Looking for sequence named '$o_entry' in blast database '$o_db'...";

my $output = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');

my @seq = qx/blastdbcmd -db $o_db -entry $o_entry/;

die "blastdbcmd could not find entry '$o_entry' in DB $o_db" if not scalar(@seq);
print STDERR " found it!\n";


my $entry_name = $seq[0];
chomp $entry_name;
$entry_name =~ s/ .*$//g;
my $seq_name = $entry_name;
my $seq = join("", @seq[1..$#seq]);
$seq =~ s/\n//g;
my ($start, $end);
if ($o_range) {
    ($start, $end) = split /-/, $o_range;
    die "problem with range definition '$o_range'" if ! $start or ! $end or $end < $start;
    $seq_name = "$seq_name:$start-$end";
}
# create BioPerl Seq object directly from scaffold sequence
my $seq_seqobj = Bio::Seq->new(-display_id => $entry_name, -seq => $seq);

$seq_seqobj = Bio::Seq->new(-display_id => $seq_name, -seq => $seq_seqobj->subseq($start, $end)) if $start and $end;

$output->write_seq($seq_seqobj);

