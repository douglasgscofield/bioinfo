#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use Statistics::Distributions;
use Data::Dumper;
use diagnostics;
#use Data::Dumper::Perltidy;
$Data::Dumper::Fill = 160;
$Data::Dumper::Indent = 1;
$Data::Dumper::Sortkeys = 1;

sub Pretty  # from http://www.perlmonks.org/?node_id=490421
{
    my @src = split(/\n/, join('', @_));
    my @dst = ();
    my $f = $Data::Dumper::Fill || 72;
    my $i = 0;
    while ($i <= $#src) {
        my $l = $src[$i];
        if (not $l =~ /[\[\{\(]\s*$/) { push(@dst, $l); $i++; next; }
        my ($p) = ($l =~ /^(\s+)/);
        my $j = $i+1;
        while ($j <= $#src) {
            my $n = $src[$j];
            my ($q) = ($n =~ /^(\s+)/);
            $n =~ s/^\s+/ /;
            if (length($l) + length($n) >= $f) { $l = $src[$i]; last; }
            $l .= $n;
            if ($q and $p and $q eq $p) { $i = $j; last; }
            $j++;
        }
        push(@dst, $l);
        $i++;
    }
    return join("\n", @dst) . "\n";
}

##use Bio::Tools::IUPAC;
# rather than depend on BioPerl, copy the hash %Bio::Tools::IUPAC::REV_IUB
# source at http://cpansearch.perl.org/src/CJFIELDS/BioPerl-1.007002/Bio/Tools/IUPAC.pm
my %REV_IUB = (
    A    => 'A',
    T    => 'T',
    U    => 'U',
    C    => 'C',
    G    => 'G',
    AC   => 'M',
    AG   => 'R',
    AT   => 'W',
    CG   => 'S',
    CT   => 'Y',
    GT   => 'K',
    ACG  => 'V',
    ACT  => 'H',
    AGT  => 'D',
    CGT  => 'B',
    ACGT => 'N',
    N    => 'N'
);

# Note: SNP Columns are transformed to 1-based, from Stacks' 0-based
# DONE: sort out >1 central snp and maxflanksnps, e.g., locus 7
# DONE: IUPAC codes for off-site SNPs, and output tag for focal SNP
# DONE: reject all templates with >2 alleles
# DONE: trim templates to 50p each side of focal SNP
# TODO: allele frequencies
# TODO: sample frequencies
# TODO: reject central SNP based on allele frequencies ?

my %LOCUS; # key, value: locus ID, stack sequence
my $n_templates = 0;

my %oland = map { $_ => 1 } qw/Gr Kv La Me Mo Na Re To Vi/;

sub process_tags_line($);
sub process_snps_line($);
sub process_sumstats_line($);
sub process_locus($);
sub condense_sumstats($);
sub evaluate_locus($);
sub print_locus_template($$);
sub dump_locus($);
sub create_locus_template($);
sub is_oland($);
sub is_central_SNP($$);

# these catalog files are produced by cstacks

my $o_dir;
my $o_snpsfile;
my $o_tagsfile;  # consensus for assembled loci
my $o_sumstatsfile;  # summary statistics per locus
my $o_alpha = 0.05;   # alpha value for frequency confidence intervals
my $o_ploidy = 2;   # ploidy of the study species
my $o_minflank = 50;
my $o_maxflanksnps = 3;
my $D_limit = 1;
my $o_templatename;
my $o_padding = 7;
my $o_verbose = 0;
my $o_hetsonly = 1; # we only want heterozygous SNPs
my $o_centralsnpgap = 5;
my $o_trimtemplate = 1;
my $o_annotatetemplate = 1;

sub usage {
    print STDERR join("", @_)."

$0: process Stacks batch files to produce SNP-assay templates
    
    --dir DIR              directory containing batch_1.{catalog.{snps,tags}.tsv.gz,sumstats.tsv} files
    --snps FILE            .snps.tsv.gz file
    --tags FILE            .tags.tsv.gz file
    --sumstats FILE        .sumstats.tsv file
    --alpha FLOAT          alpha (two-sided) for SNP frequency confidence intervals via Wilson's method [$o_alpha]
    --minflank INT         minimum flank length; central SNPs are those between flanks [$o_minflank]
    --maxflanksnps INT     maximum number of het SNPs in each flank [$o_maxflanksnps]
    --centralsnpgap INT    minimum gap (bp) between central SNPs [$o_centralsnpgap]
    --templatename STRING  name prefix for templates
    --padding INT          0-padded with for number suffix on template name [$o_padding]
    --trimtemplate         trim templates to have minflank+SNP+minflank bases
    --annotatetemplate     add annotation columns to template output, for filtering
    --D_limit INT          only read this number of templates; for debugging
    --verbose              verbose debugging output

";
    exit 1;
}

usage() if ! scalar(@ARGV);

GetOptions("dir=s" => \$o_dir,
           "snps=s" => \$o_snpsfile,
           "tags=s" => \$o_tagsfile,
           "sumstats=s" => \$o_sumstatsfile,
           "minflank=i" => \$o_minflank,
           "maxflanksnps=i" => \$o_maxflanksnps,
           "centralsnpgap=i" => \$o_centralsnpgap,
           "templatename=s" => \$o_templatename,
           "padding=i" => \$o_padding,
           "trimtemplate" => \$o_trimtemplate,
           "annotatetemplate" => \$o_annotatetemplate,
           "D_limit=i" => \$D_limit,
           "verbose:1" => \$o_verbose,
) or usage();

if ($o_dir) {
    $o_snpsfile ||= "$o_dir/batch_1.catalog.snps.tsv.gz";
    $o_tagsfile ||= "$o_dir/batch_1.catalog.tags.tsv.gz";
    $o_sumstatsfile ||= "$o_dir/batch_1.sumstats.tsv";
    $o_templatename ||= "${o_dir}_";
} else {
    $o_snpsfile ||= "batch_1.catalog.snps.tsv.gz";
    $o_tagsfile ||= "batch_1.catalog.tags.tsv.gz";
    $o_sumstatsfile ||= "batch_1.sumstats.tsv";
    $o_templatename ||= "batch_1_";
}

open my $f_snpsfile, '-|', "gzip -dcf $o_snpsfile" or die "could not open $o_snpsfile for reading: $!";
unless ((my $nop = <$f_snpsfile>) =~ /^# cstacks/) { die "First line of snps file $o_snpsfile not '# cstacks ...'"; }
#die "First line of snps file $o_snpsfile not '# cstacks ...'" unless (my $nop = <$f_snpsfile>) =~ /^# cstacks/;

open my $f_tagsfile, '-|', "gzip -dcf $o_tagsfile" or die "could not open $o_tagsfile for reading: $!";
unless ((my $nop = <$f_tagsfile>) =~ /^# cstacks/) { die "First line of tags file $o_tagsfile not '# cstacks ...'"; }

open my $f_sumstatsfile, '-|', "gzip -dcf $o_sumstatsfile" or die "could not open $o_sumstatsfile for reading: $!";
unless ((my $nop = <$f_sumstatsfile>) =~ /^# 1\t/) { die "First line of tags file $o_sumstatsfile not '# 1<tab>...'"; }


# Read the tags first to fill in %LOCUS, one tag per consensus stack
while (<$f_tagsfile>) {
    my $h = process_tags_line($_);
    $LOCUS{$h->{locus}} = $h;
    last if $D_limit && $. > $D_limit;
}


# Read SNPs, add them to %LOCUS
while (<$f_snpsfile>) {  # snps, 0, 1, or more than 1 per consensus stack
    my $h = process_snps_line($_);
    last if ! exists $LOCUS{$h->{locus}};  # if we have not read this locus, stop
    push @{$LOCUS{$h->{locus}}->{snps}}, $h;
}
#say Dumper(\%LOCUS);

# Read sumstats information for each SNP, add to locus if central SNP
while (<$f_sumstatsfile>) {  # summary statistics for each SNP in each population
    next if /^#/;
    my $h = process_sumstats_line($_);
    #print "$o_sumstatsfile:$.: loc:$h->{locus} \@$h->{column} $h->{pop}:$h->{nindiv} $h->{alp}/$h->{alq} ".sprintf("%.3f", $h->{freqp}).($h->{oland}?" oland":"").($h->{private}?" private":"") if 1;
    my $out = "$o_sumstatsfile:$.: loc:$h->{locus} \@$h->{column} $h->{pop}:$h->{nindiv} $h->{alp}/$h->{alq} ".sprintf("%.3f", $h->{freqp});
    if (! exists $LOCUS{$h->{locus}}) { # if we have not read this locus, stop
        say "$out LOCUS NOT READ, BAILING";
        last;
    }
    if (! is_central_SNP($LOCUS{$h->{locus}}, $h->{column})) {
        #print "\n" if 1;
        next;
    }
    say "$out CENTRAL" if 1;
    push @{$LOCUS{$h->{locus}}->{sumstats}->{$h->{column}}}, $h;
}
say Pretty(Dumper(\%LOCUS));
#exit 1;


process_locus($LOCUS{$_}) foreach sort { $a <=> $b } keys %LOCUS;

if ($o_annotatetemplate) {
    say STDOUT "*** Output: name locus focal_SNP_0-based n_lsnps:n_rsnps template";
} else {
    say STDOUT "*** Output: name template";
}
say STDOUT "*** Templates for $n_templates loci".($o_trimtemplate ? " (flanks trimmed to $o_minflank bp)" : "");


###############################
#
# subroutines
#

sub process_locus($) {
    my $locus = shift;
    condense_sumstats($locus);
    if (evaluate_locus($locus)) {
        create_locus_template($locus);
        print_locus_template($locus, $o_annotatetemplate);
    }
    dump_locus($locus) if $o_verbose;
}


sub process_tags_line($) {
    chomp $_[0];
    my (undef, $sample, $locus, undef, undef, undef, $type, $component, $seqid, $seq, undef, undef, undef, undef) = split(/\t/, $_[0]);
    die "$o_tagsfile:$.: sequence type '$type' not 'consensus'" if $type ne 'consensus';
    #return { locus => $locus, seq => $seq, seqid => $seqid, seqlen => length($seq) };
    return { locus => $locus, seq => $seq, seqlen => length($seq) };
}


sub process_snps_line($) {
    chomp $_[0];
    my (undef, $sample, $locus, $column, $type, $LR, $al0, $al1, $al2, $al3) = split(/\t/, $_[0]);
    die "$o_snpsfile:$.: snp type '$type' not 'E'" if $type ne 'E';
    $column += 1;  # convert to 1-based
    return { locus => $locus, column => $column, type => $type, al0 => $al0, al1 => $al1, al2 => $al2, al3 => $al3 };
}


sub process_sumstats_line($) {
    chomp $_[0];
    #   $batch, $locus, $chrom, $bp_cum, $column, $pop, $alp, $alq, $nindiv, $pfreq, $obshet, $obshom, $exphet, $exphom, $pi,   $pismooth, $P_pismooth, $Fis , $Fis_smooth, $P_Fis_smooth, $private
    my ($batch, $locus, $chrom, undef  , $column, $pop, $alp, $alq, $nindiv, $freqp, undef  , undef  , undef  , undef  , undef, undef    , undef      , undef, undef      , undef        , undef) = split(/\t/, $_[0]);
    die "$o_sumstatsfile:$.: Chr '$chrom' not 'un'" if $chrom ne 'un';
    $column += 1;  # convert to 1-based
    return { locus => $locus, column => $column, pop => $pop, nindiv => $nindiv, alp => $alp, alq => $alq, freqp => $freqp };
}


sub is_oland($) {
    return exists($oland{substr($_[0], 0, 2)});
}


sub is_central_SNP($$) {
    my ($locus, $column) = @_;
    return ($column > $o_minflank and $column < ($locus->{seqlen} - $o_minflank + 1)) ? 1 : 0;
}

sub condense_sumstats($) {
    my $locus = shift;
    my $sumstats = $locus->{sumstats};
    foreach my $k (sort keys %$sumstats) {
        say "locus $locus->{locus} \@ $k";
        my $npop = @{$sumstats->{$k}};
        my %oland;
        my %other;
        foreach (@{$sumstats->{$k}}) {
            my $nindiv  = $_->{nindiv};
            my $nhap    = $nindiv * $o_ploidy;
            my $nq      = $nhap - ($nhap * $_->{freqp});
            my $is_poly = $_->{freqp} < 1.0 && $_->{freqp} > 0.0;
            my $update_group = sub {
                my $h = shift;
                $h->{npop}++;
                $h->{nindiv}  += $nindiv;
                $h->{nhap}    += $nhap;
                $h->{nq}      += $nq;
                $h->{npoly}   += $is_poly;
                $h->{qfreq} = $h->{nq} / $h->{nhap};
                $h->{qci} = [ confidence_interval($h->{nhap}, $h->{nq}, $o_alpha, 0) ];
            };
            $update_group->(is_oland($_->{pop}) ? \%oland : \%other);
        }
        print "oland ";
        say "npop $_->{npop} nindiv $_->{nindiv} nhap $_->{nhap} nq $_->{nq} npoly $_->{npoly}\t$_->{qfreq}\t[$_->{qci}->[0], $_->{qci}->[1]]" for (\%oland);
        print "other ";
        say "npop $_->{npop} nindiv $_->{nindiv} nhap $_->{nhap} nq $_->{nq} npoly $_->{npoly}\t$_->{qfreq}\t[$_->{qci}->[0], $_->{qci}->[1]]" for (\%other);
    }
    # count oland and non-oland populations that are polymorphic
    # calculate frequencies
    # calculate binomial probabilities of observed frequencies, with confidence intervals
    # compute confidence intervals using Wilson
}

sub confidence_interval($$$$) {
    # https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    my ($n, $success, $alpha, $method) = @_;
    my $z = Statistics::Distributions::udistr($alpha / 2);  # two-sized z-score
    usage("alpha $alpha seems inappropriate, z-score is $z") if $z < 0;
    my $p = $success / $n;
    my @ci;
    if (! $method) { # wilson's method
        my $z2 = $z * $z;
        my $pm = $z * sqrt(((1.0/$n) * $success * ($n - $success)) + 0.25*$z2);
        @ci = map { (1.0/($n + $z2)) * ($success + 0.5*$z2 + ($_ * $pm)) } (-1.0, +1.0);
    } elsif ($method == 1) { # adjust p and use normal approximation
        $p = ($success + 2) / ($n + 4);
        my $pm = $z * sqrt((1.0 / $n) * $p * (1.0 - $p));  # the +/- part
        @ci = map { $p + ($_ * $pm) } (-1.0, +1.0);
    } elsif ($method == 2) { # normal approximation, alternative factorisation
        my $pm = $z * sqrt((1.0 / $n) * $success * ($n - $success));  # the +/- part
        @ci = map { (1.0 / $n) * ($success + ($_ * $pm)) } (-1.0, +1.0);
    } else { die "unknown method: $method"; }
    @ci = map { sprintf("%.4f", $_) } @ci;
    return @ci;
}

sub print_locus_template($$) {
    my ($locus, $do_annotate) = @_;
    my $name = $locus->{name};
    my $template = $locus->{template};
    if ($o_trimtemplate) {
        my $x;
        $x = $locus->{focal_SNP} - 1 - $o_minflank;
        substr($template, 0, $x) = "" if $x;
        $x = $locus->{seqlen} - $locus->{focal_SNP} - $o_minflank;
        substr($template, -$x) = "" if $x;
    }
    if ($do_annotate) {
        say STDOUT join("\t", $name, $locus->{locus}, $locus->{focal_SNP} - 1, "$locus->{n_lsnps}:$locus->{n_rsnps}", $template);
    } else {
        say STDOUT "$name $template";
    }
    ++$n_templates;
}

sub dump_locus($) {
    my $locus = shift;
    say STDOUT "locus $locus->{locus} consensus $locus->{seqlen} bp:";
    say STDOUT "$locus->{seq}";
    if (exists $locus->{snps}) {
        my $focal = exists $locus->{focal_SNP} ?  $locus->{focal_SNP} : 0;
        say STDOUT (' ' x ($_->{column} - 1)).($_->{column} == $focal ? "*" : "^")."$_->{al0}$_->{al1}$_->{al2}".'@'."$_->{column}" foreach @{$locus->{snps}};
    } else {
        say STDOUT "-- no snps";
    }
    if (exists $locus->{template}) {
        print_locus_template($locus, $o_annotatetemplate);
    } else {
        say STDOUT "-- no template";
    }
}


sub create_locus_template($) {
    my $locus = shift;
    if (! exists $locus->{focal_SNP}) {
        say "no focal SNP for locus $locus->{locus}, not creating template" if $o_verbose;
        return;
    }
    if ($locus->{focal_SNP} == 0) {
        say "*** focal_SNP key exists but value == 0";
        exit 1;
    }
    my $name = $o_templatename.sprintf("%.${o_padding}d", $locus->{locus});
    my $template = $locus->{seq};
    my $focal_SNP;
    # first, non-focal SNPs get converted to IUPAC symbols; these substitutions do not change template length
    foreach (@{$locus->{snps}}) {
        if ($_->{column} == $locus->{focal_SNP}) {
            $focal_SNP = $_;
            next;
        }
        substr($template, $_->{column} - 1, 1) = $REV_IUB{uc join '', sort ( $_->{al0}, $_->{al1} )};
    }
    # now the focal SNP; changes template length
    my @a = sort map { uc } ( $focal_SNP->{al0}, $focal_SNP->{al1} );
    substr($template, $focal_SNP->{column} - 1, 1) = "[$a[0]/$a[1]]";
    # pack up and go
    $locus->{name} = $name;
    $locus->{template} = $template;
}


sub evaluate_locus($) {
    # is this a locus suitable for creating a probe?
    my $locus = shift;
    print STDOUT "evaluating locus $locus->{locus} ... " if $o_verbose;
    if (! exists $locus->{snps}) {
        say STDOUT "no SNPs" if $o_verbose; 
        return 0;
    }

    # sort SNPs by column
    $locus->{snps} = [ sort { $a->{column} <=> $b->{column} } @{$locus->{snps}} ];
    ##my @s = sort { $a->{column} <=> $b->{column} } @{$LOCUS{$l}->{snps}};

    # reject templates with any 3+ allele SNPs
    my @snps3 = grep { $_->{al3} ne '-' } @{$locus->{snps}};
    if (@snps3) {
        say STDOUT 'at least one 3+ allele SNP' if $o_verbose; 
        return 0;
    }

    # reject template with no central SNP
    #my @csnps = grep { $_->{column} > $o_minflank and $_->{column} < ($locus->{seqlen} - $o_minflank + 1) } @{$locus->{snps}};
    my @csnps = grep { is_central_SNP($locus, $_->{column}) } @{$locus->{snps}};
    if (! @csnps) {
        say STDOUT 'no central snps' if $o_verbose;
        return 0;
    }

    # get SNPs on the flanks
    my @lsnps = grep { $_->{column} <= $o_minflank } @{$locus->{snps}};
    my $l_column = @lsnps > 0 ? $lsnps[$#lsnps]->{column} : 0;  # column of the rightmost left SNP
    my @rsnps = grep { $_->{column} >= ($locus->{seqlen} - $o_minflank + 1) } @{$locus->{snps}};
    my $r_column = @rsnps > 0 ? $rsnps[0]->{column} : 0;  # column of the leftmost right SNP

    if (@csnps == 1) {
        if (($l_column && ($csnps[0]->{column} - $l_column) < $o_centralsnpgap) or
            ($r_column && ($r_column - $csnps[0]->{column}) < $o_centralsnpgap)) {
            say STDOUT 'single central SNP too close to SNP in flank' if $o_verbose;
            return 0;
        }
    } else {  # more than one central SNP
        my @close = grep { ($csnps[$_ + 1]->{column} - $csnps[$_]->{column}) < $o_centralsnpgap }  0 .. ($#csnps - 1);
        if (@close) {
            say STDOUT 'at least two central SNPs are too close to each other' if $o_verbose;
            return 0;
        }
        if (($l_column && ($csnps[0]->{column} - $l_column) < $o_centralsnpgap) or
            ($r_column && ($r_column - $csnps[$#csnps]->{column}) < $o_centralsnpgap)) {
            say STDOUT 'at least one central SNP too close to SNP in flank' if $o_verbose;
            return 0;
        }
        # choose a central SNP, move the other(s) to the flanks
        # 1) favour equal numbers of SNPs on both sides of central SNP
        # 2) favour a more central SNP for the focal SNP over a less central SNP
        while (@csnps > 1) {
            if (@lsnps < @rsnps) { # move leftmost central SNP to end of left flank
                push @lsnps, shift @csnps;
            } elsif (@lsnps > @rsnps) { # move rightmost central SNP to beginning of right flank
                unshift @rsnps, pop @csnps;
            } else {
                # we have the same number of SNPs on the flanks, either 0 or >0
                # keep the 'most' central of the two outer SNPs, or a random one if both are equally central
                my $lbias = abs(1.0 - ($csnps[0]->{column} / ($locus->{seqlen} - $csnps[0]->{column} + 1)));
                my $rbias = abs(1.0 - ($csnps[$#csnps]->{column} / ($locus->{seqlen} - $csnps[$#csnps]->{column} + 1)));
                if ($lbias > $rbias or ($lbias == $rbias && rand() < 0.5)) { # left is less central OR equally so and random draw
                    push @lsnps, shift @csnps;
                } else {  # right is less central OR equally so and failed random draw
                    unshift @rsnps, pop @csnps;
                }
            }
        }
    }
    if (@lsnps > $o_maxflanksnps) {
        say STDOUT 'too many SNPs on left flank' if $o_verbose;
        return 0;
    }
    if (@rsnps > $o_maxflanksnps) {
        say STDOUT 'too many SNPs on right flank' if $o_verbose;
        return 0;
    }

    say STDOUT "SUITABLE SNP at column $csnps[0]->{column}" if $o_verbose;

    $locus->{n_lsnps} = scalar @lsnps;
    $locus->{n_rsnps} = scalar @rsnps;
    $locus->{focal_SNP} = $csnps[0]->{column};
    return $locus->{focal_SNP};
}


__END__

foreach (keys %oland) {
    say;
    say "$_ ".(is_oland($_) ? 'Oland' : 'nope');
}
