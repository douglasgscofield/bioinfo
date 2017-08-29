#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use Data::Dumper::Perltidy;

# Note: SNP Columns are transformed to 1-based, from Stacks' 0-based
# TODO: sort out >1 central snp and maxflanksnps, e.g., locus 7

my %LOCUS; # key, value: locus ID, stack sequence

my %oland = map { $_ => 1 } qw/Gr Kv La Me Mo Na Re To Vi/;

sub is_oland($) {
    return exists($oland{substr($_[0], 0, 2)});
}

# make sure SNPs file is good

# these catalog files are produced by cstacks

my $o_snpsfile = "batch_1.catalog.snps.tsv.gz";
my $o_tagsfile = "batch_1.catalog.tags.tsv.gz";  # consensus for assembled loci
my $o_minflank = 50;
my $o_maxflanksnps = 2;
my $D_limit = 99;
my $o_hetsonly = 1; # we only want heterozygous SNPs

GetOptions("snps=s" => \$o_snpsfile,
           "tags=s" => \$o_tagsfile,
           "minflank=i" => \$o_minflank,
           "maxflanksnps=i" => \$o_maxflanksnps,
           "D_limit=i" => \$D_limit,
) or die "unrecognized option: $!";


open my $f_snpsfile, '-|', "gzip -dcf $o_snpsfile" or die "could not open $o_snpsfile for reading: $!";
unless ((my $nop = <$f_snpsfile>) =~ /^# cstacks/) { die "First line of snps file $o_snpsfile not '# cstacks ...'"; }
#die "First line of snps file $o_snpsfile not '# cstacks ...'" unless (my $nop = <$f_snpsfile>) =~ /^# cstacks/;

open my $f_tagsfile, '-|', "gzip -dcf $o_tagsfile" or die "could not open $o_tagsfile for reading: $!";
unless ((my $nop = <$f_tagsfile>) =~ /^# cstacks/) { die "First line of tags file $o_tagsfile not '# cstacks ...'"; }

sub process_tags_line($) {
    chomp $_[0];
    my (undef, $sample, $locus, undef, undef, undef, $type, $component, $seqid, $seq, undef, undef, undef, undef) = split(/\t/, $_[0]);
    die "$o_tagsfile:$.: sequence type '$type' not 'consensus'" if $type ne 'consensus';
    return { locus => $locus, seq => $seq, seqid => $seqid, seqlen => length($seq) };
}
sub process_snps_line($) {
    chomp $_[0];
    my (undef, $sample, $locus, $column, $type, $LR, $al0, $al1, $al2, $al3) = split(/\t/, $_[0]);
    die "$o_snpsfile:$.: snp type '$type' not 'E'" if $type ne 'E';
    return { locus => $locus, column => $column + 1, type => $type, al0 => $al0, al1 => $al1, al2 => $al2, al3 => $al3 };
}
sub dump_locus($) {
    my $l = shift;
    say STDOUT "locus $l consensus $LOCUS{$l}->{seqlen} bp:\n$LOCUS{$l}->{seq}";
    if (exists $LOCUS{$l}->{snps}) {
        my @s = sort { $a->{column} <=> $b->{column} } @{$LOCUS{$l}->{snps}};
        say STDOUT (' ' x ($_->{column} - 1))."^$_->{al0}$_->{al1}$_->{al2} c$_->{column}" foreach @s;
        evaluate_locus($l);
    } else {
        say "-- no snps";
    }
}
sub evaluate_locus($) {
    my $l = shift; my $seqlen = $LOCUS{$l}->{seqlen};
    # is this a locus suitable for creating a probe?  all SNPs must be two-allele SNPs
    print STDOUT "evaluating locus $l ... ";
    my @csnps = grep { $_->{column} > $o_minflank and $_->{column} < ($seqlen - $o_minflank + 1) and $_->{al3} eq '-' } @{$LOCUS{$l}->{snps}};
    if (! @csnps) { say STDOUT 'no central snps'; return; }
    my @lsnps = grep { $_->{column} <= $o_minflank and $_->{al3} eq '-' } @{$LOCUS{$l}->{snps}};
    if (@lsnps > $o_maxflanksnps) { say STDOUT 'too many snps on left flank'; return; }
    my @rsnps = grep { $_->{column} >= ($seqlen - $o_minflank + 1) and $_->{al3} eq '-' } @{$LOCUS{$l}->{snps}};
    if (@rsnps > $o_maxflanksnps) { say STDOUT 'too many snps on right flank'; return; }
    say STDOUT 'SUITABLE';
}
while (<$f_tagsfile>) {  # tags first, one tag per consensus stack
    my $h = process_tags_line($_);
    $LOCUS{$h->{locus}} = $h;
    last if $D_limit && $. > $D_limit;
}
while (<$f_snpsfile>) {  # snps, 0, 1, or more than 1 per consensus stack
    my $h = process_snps_line($_);
    last if ! exists $LOCUS{$h->{locus}};  # if we have not read this locus, stop
    push @{$LOCUS{$h->{locus}}->{snps}}, $h;
}
#say Dumper(\%LOCUS);

dump_locus($_) foreach sort { $a <=> $b } keys %LOCUS;


__END__

foreach (keys %oland) {
    say;
    say "$_ ".(is_oland($_) ? 'Oland' : 'nope');
}
