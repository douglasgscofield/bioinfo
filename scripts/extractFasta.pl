#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Seq;  # the BioPerl sequence manipulation package
use Bio::SeqIO;  # the BioPerl sequence file manipulation package
use Getopt::Long;

my $o_db;
my $o_entry;
my $o_range;
my $o_strict = 0;
my $endpoints_swapped = 0;
my $o_help;

my $USAGE = "
Extract a sequence or subsequence from a Blast database

USAGE: $0 --db database --entry fasta-sequence-name [ --range L-R ]

First build the blast database using

    makeblastdb -parse_seqids -in sequence.fa -dbtype nucl|prot

with -dbtype dependent upon your input data.

Then you can extract the (sub)sequence of interest.  This is printed
to stdout in Fasta format.

If a subsequence is requested with the --range option, the range
of the subsequence is added to the sequence identifier in the output.

OPTIONS:

    --db blast-db       The same name you would give to blastdbcmd,
                        which is exactly what this script does

    --entry sequence    Name of the sequence to extract from blast-db

    --range LOW-HIGH    1-based positions of a subsequence to extract
                        from the sequence, once it is found.  If HIGH < LOW
                        the ends of the range are swapped unless --strict

    --strict            do not be tolerant of HIGH < LOW in --range

";

GetOptions("db=s"    => \$o_db,
           "entry=s" => \$o_entry,
           "range=s" => \$o_range,
           "strict"  => \$o_strict,
           "help"    => \$o_help) or die $USAGE;
if ($o_help or not $o_db or not $o_entry) { print STDERR $USAGE; exit 0; }

my ($start, $end);
if ($o_range) {
    ($start, $end) = split /-/, $o_range;
    die "either end of range is 0 or unknown '$o_range'" if ! $start or ! $end;
    if ($end < $start and not $o_strict) {
        ($start, $end) = ($end, $start);
        $endpoints_swapped = 1;
    }
    die "order of range extents incorrect '$o_range'" if $end < $start;
}

print STDERR "Looking for sequence named '$o_entry' in blast database '$o_db'...";

my $output = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');

my @seq = qx/blastdbcmd -db $o_db -entry '$o_entry'/;

die "blastdbcmd could not find entry '$o_entry' in DB $o_db" if not scalar(@seq);
print STDERR " found it!" . ($endpoints_swapped ? "  range endpoints swapped" : "") . "\n";


my $entry_name = $seq[0];
chomp $entry_name;
$entry_name =~ s/^>//g;  # remove starting >
$entry_name =~ s/ .*$//g;  # remove everything after the first space, keeping only the identifier
my $seq_name = $entry_name;
my $seq = join("", @seq[1..$#seq]);
$seq =~ s/\n//g;
if ($o_range) {
    $seq_name = "$seq_name:$start-$end";
}
# create BioPerl Seq object directly from scaffold sequence
my $seq_seqobj = Bio::Seq->new(-display_id => $entry_name, -seq => $seq);

$seq_seqobj = Bio::Seq->new(-display_id => $seq_name, -seq => $seq_seqobj->subseq($start, $end)) if $start and $end;

$output->write_seq($seq_seqobj);

