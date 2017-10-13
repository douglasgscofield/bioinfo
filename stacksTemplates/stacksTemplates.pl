#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use Statistics::Distributions;
use Data::Dumper;
#use diagnostics;
#use Data::Dumper::Perltidy;
$Data::Dumper::Fill = 160;
$Data::Dumper::Indent = 1;
$Data::Dumper::Sortkeys = 1;

sub Pretty  # condense Dumper output, from http://www.perlmonks.org/?node_id=490421
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
    A    => 'A', T    => 'T', U    => 'U', C    => 'C', G    => 'G',
    AC   => 'M', AG   => 'R', AT   => 'W', CG   => 'S', CT   => 'Y', GT   => 'K',
    ACG  => 'V', ACT  => 'H', AGT  => 'D', CGT  => 'B', ACGT => 'N', N    => 'N'
);

# Note: SNP Columns are transformed to 1-based, from Stacks' 0-based

my %LOCUS;
my $n_templates = 0;

my %NORTH = map { $_ => 1 } qw/Ab Hs Li Ru Sk Tr/;
my %OLAND = map { $_ => 1 } qw/Gr Kv La Me Mo Na Re To Vi/;
my %REGIONS = map { $_ => $_ } qw/north other oland/;
my %POPULATIONS_SEEN;

sub is_oland($) {
    return exists($OLAND{substr($_[0], 0, 2)});
}

sub pop_region($) {
    my ($pop, $region) = (substr($_[0], 0, 2), "");
    if (exists($OLAND{$pop})) {
        $region = $REGIONS{oland};
    } elsif (exists($NORTH{$pop})) {
        $region = $REGIONS{north};
    } else {
        $region = $REGIONS{other};
    }
    return $region;
}


sub process_tags_line($);
sub process_snps_line($);
sub process_sumstats_line($);
sub process_locus($);
sub condense_sumstats($);
sub evaluate_locus($);
sub print_locus_template($$);
sub dump_locus($);
sub create_locus_template($);
sub is_central_SNP($$);
sub binomial_confint($$$$);
sub check_SNP_has_sumstats($$);
sub check_SNP_sumstats_OK($$);

# these catalog files are produced by cstacks

my $o_dir;
my $o_snpsfile;
my $o_tagsfile;  # consensus for assembled loci
my $o_sumstatsfile;  # summary statistics per locus
my $o_alpha = 0.05;   # alpha value for frequency confidence intervals
my $o_ploidy = 2;   # ploidy of the study species
my $o_minflank = 50;
my $o_maxflanksnps = 3;
my $D_limit = 0;
my $o_name;
my $o_extended;
my $o_padding = 7;
my $o_verbose = 0;
my $o_centralsnpgap = 5;
my $o_trimtemplate = 1;
my $o_annotatetemplate = 1;
my $o_freqdigits = 4;

# criteria of central SNP is OK to use
my $o_crit_n_region = 2;
my $o_crit_nwithq = 1;
my $o_crit_qfreq = 0.01;

my $N_locus_no_SNP = 0;
my $N_locus_SNP_3allele = 0;
my $N_locus_no_central_SNP = 0;
my $N_locus_central_SNP_no_sumstats = 0;
my $N_locus_central_SNP_no_sumstats_OK = 0;
my $N_locus_central_SNP_crowded = 0;
my $N_locus_central_SNP_multiple_one_OK = 0;
my $N_locus_central_SNP_multiple_OK = 0;
my $N_locus_central_SNP_multiple_remaining = 0;
my $N_locus_flank_SNP_toomany = 0;

sub usage {
    print STDERR join("", @_)."

$0: process Stacks batch files to produce SNP-assay templates

Input files are given with these options.  The --dir option can be used alone.
    
    --dir DIR              directory containing batch_1.{catalog.{snps,tags}.tsv.gz,sumstats.tsv} files
    --snps FILE            .snps.tsv.gz file
    --tags FILE            .tags.tsv.gz file
    --sumstats FILE        .sumstats.tsv file

To be considered, stack consensus sequences must pass these basic criteria:

    --minflank INT         minimum flank length; central SNPs are those between flanks [$o_minflank]
    --maxflanksnps INT     maximum number of SNPs in each flank [$o_maxflanksnps]
    --centralsnpgap INT    minimum gap (bp) between central SNPs [$o_centralsnpgap]

Only central SNPs that pass simple filtering according to the following criteria are considered:

    --crit_n_region INT    min number of regions stack must be scored in [$o_crit_n_region]
    --crit_nwithq INT      min number of populations with the minor allele [$o_crit_nwithq]
    --crit_qfreq FLOAT     min minor allele frequency across all scored populations [$o_crit_qfreq]

For all SNP minor alleles, binomial confints are calculated as well:

    --alpha FLOAT          alpha for minor allele frequency confint (Wilson's method) [$o_alpha]

All central SNPs that are considered are then ranked using these criteria, in order:

    Highest number of populations with minor allele across all regions
    Highest number of regions with minor allele
    Highest number of regions stack was scored in
    Highest lower bound of minor allele frequency confint across all regions
    Highest minor allele frequency across all regions

The central SNP that is ranked most highly is used, and the others are shifted to the flanks.

These options control output.

    --name STRING          name prefix for templates, only one of --name or --extended allowed
    --extended STRING      use extended template names, with STRING for the nonvariable prefix
                           Extended names are particular to this study.
                           STRINGR#cP##p###_locusnumber0padded
                           STRING: the prefix
                           R#c   : the regions containing the minor allele
                                   R3A : all regions
                                   R1N : only north
                                   R1T : only other
                                   R1L : only oland
                                   R2X : north and other
                                   R2Y : north and oland
                                   R2Z : other and oland
                           P##   : Total number of populations with minor allele, 0-padded
                           p###  : Number north, other, oland populations containing minor allele
    --padding INT          0-padded with for number suffix on template name [$o_padding]
    --trimtemplate         trim templates to have minflank+SNP+minflank bases [$o_trimtemplate]
    --annotatetemplate     add annotation columns to template output, for filtering [$o_annotatetemplate]
    --freqdigits INT       digits right of decimal point for allele frequences [$o_freqdigits]
    --D_limit INT          only read this number of templates; for debugging
    --verbose              verbose debugging output

";
    exit 1;
}

usage() if ! scalar(@ARGV);

GetOptions("dir=s"             => \$o_dir,
           "snps=s"            => \$o_snpsfile,
           "tags=s"            => \$o_tagsfile,
           "sumstats=s"        => \$o_sumstatsfile,
           "minflank=i"        => \$o_minflank,
           "maxflanksnps=i"    => \$o_maxflanksnps,
           "centralsnpgap=i"   => \$o_centralsnpgap,
           "crit_n_region=i"   => \$o_crit_n_region,
           "crit_nwithq=i"     => \$o_crit_nwithq,
           "crit_qfreq=f"      => \$o_crit_qfreq,
           "alpha=f"           => \$o_alpha,
           "name=s"            => \$o_name,
           "extended=s"        => \$o_extended,
           "padding=i"         => \$o_padding,
           "trimtemplate"      => \$o_trimtemplate,
           "annotatetemplate"  => \$o_annotatetemplate,
           "freqdigits=i"      => \$o_freqdigits,
           "D_limit=i"         => \$D_limit,
           "verbose:1"         => \$o_verbose,
) or usage();

($o_name xor $o_extended) or die "only one of --name and --extended allowed";

if ($o_dir) {
    $o_snpsfile ||= "$o_dir/batch_1.catalog.snps.tsv.gz";
    $o_tagsfile ||= "$o_dir/batch_1.catalog.tags.tsv.gz";
    $o_sumstatsfile ||= "$o_dir/batch_1.sumstats.tsv";
    $o_name ||= "${o_dir}_";
} else {
    $o_snpsfile ||= "batch_1.catalog.snps.tsv.gz";
    $o_tagsfile ||= "batch_1.catalog.tags.tsv.gz";
    $o_sumstatsfile ||= "batch_1.sumstats.tsv";
    $o_name ||= "batch_1_";
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
    $POPULATIONS_SEEN{$h->{pop}}++;
    my $out = "$o_sumstatsfile:$.: loc:$h->{locus} \@$h->{column} $h->{pop}:$h->{nindiv} $h->{alp}/$h->{alq} ".sprintf("%.${o_freqdigits}f", $h->{freqp});
    if (! exists $LOCUS{$h->{locus}}) { # if we have not read this locus, stop
        say "$out LOCUS NOT READ, SKIPPING" if $o_verbose;
        next;
    }
    next if ! is_central_SNP($LOCUS{$h->{locus}}, $h->{column});
    say "$out CENTRAL" if $o_verbose;
    push @{$LOCUS{$h->{locus}}->{sumstats}->{$h->{column}}}, $h;
}
say Pretty(Dumper(\%LOCUS)) if $o_verbose;
#exit 1;


process_locus($LOCUS{$_}) foreach sort { $a <=> $b } keys %LOCUS;

say STDERR "
*** dir . . . . . . . . . . $o_dir
*** SNPs file               $o_snpsfile
*** tags file . . . . . . . $o_tagsfile
*** sumstats file           $o_sumstatsfile
*** template name . . . . . $o_name
*** extended template name  $o_extended
*** frequency digits. . . . $o_freqdigits
*** locus number padding    $o_padding
*** alpha for confint . . . $o_alpha
*** ploidy                  $o_ploidy
*** D_limit (debugging)   . $D_limit
*** verbose (debugging) . . $o_verbose
***
*** minimum flank length                             $o_minflank
*** maximum SNPs on each flank . . . . . . . . . . . $o_maxflanksnps
*** minimum gap between central SNPs                 $o_centralsnpgap
*** minimum number of scored regions . . . . . . . . $o_crit_n_region
*** minimum number of populations with minor allele  $o_crit_nwithq
*** minimum minor allele frequency . . . . . . . . . $o_crit_qfreq
***
*** Templates: ".($o_extended ? "extended_name" : "name")." ".($o_annotatetemplate ?  "regions_scored(pops nor,oth,ola) regions_withq(pops nor,oth,ola) qfreq(nor,oth,ola) SNP_column n_lsnps:n_rsnps" : "")." template".
($o_extended ? "
***
*** extended_name is <prefix>R#cP##p###_locusnumber0padded:
***     <prefix> the argument to --extended
***     R#c is the regions the minor allele was seen in, with c the region code
***         R1N: north, R1T: other, R1L: oland
***         R2X: north and other, R2Y: north and oland, R2Z: other and oland
***         R3A: all three
***     P## is the total number of populations with the minor allele, 0-padded
***     p### is the number of populations the minor allele was seen in
***      ^   number of populations in north (single digit)
***       ^  number of populations in other (single digit)
***        ^ number of populations in oland (single digit)" : "")."
***
*** Templates for $n_templates loci".($o_trimtemplate ? " (flanks trimmed to $o_minflank bp)" : "")."
***
*** N_locus_no_SNP                          $N_locus_no_SNP
*** N_locus_SNP_3allele . . . . . . . . . . $N_locus_SNP_3allele
*** N_locus_no_central_SNP                  $N_locus_no_central_SNP
*** N_locus_central_SNP_no_sumstats . . . . $N_locus_central_SNP_no_sumstats
*** N_locus_central_SNP_no_sumstats_OK      $N_locus_central_SNP_no_sumstats_OK
*** N_locus_central_SNP_crowded . . . . . . $N_locus_central_SNP_crowded
*** N_locus_central_SNP_multiple_one_OK     $N_locus_central_SNP_multiple_one_OK
*** N_locus_central_SNP_multiple_OK . . . . $N_locus_central_SNP_multiple_OK
*** N_locus_central_SNP_multiple_remaining  $N_locus_central_SNP_multiple_remaining
*** N_locus_flank_SNP_toomany . . . . . . . $N_locus_flank_SNP_toomany
***
*** The following populations were observed (name, region, number of scored stacks)";
foreach (sort keys %POPULATIONS_SEEN) {
    say STDERR "*** ".join("\t", $_, pop_region($_), $POPULATIONS_SEEN{$_});
}


###############################
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


sub is_central_SNP($$) {
    my ($locus, $column) = @_;
    return ($column > $o_minflank and $column < ($locus->{seqlen} - $o_minflank + 1)) ? 1 : 0;
}


sub condense_sumstats($) {
    my $locus = shift;
    my $sumstats = $locus->{sumstats};
    foreach my $k (sort keys %$sumstats) {
        say "condense_sumstats: locus $locus->{locus} \@ $k" if $o_verbose;
        my $ans = { };
        # initialize totals
        $ans->{_npop}   = scalar(@{$sumstats->{$k}});
        $ans->{_nindiv} = $ans->{_nhap} = $ans->{_nq} = $ans->{_nwithq} = $ans->{_qfreq} = 0;
        $ans->{_qci}    = [ (0, 0) ];
        foreach (keys %REGIONS) {  # initialise for all regions
            $ans->{$_}->{npop} = $ans->{$_}->{nindiv} = $ans->{$_}->{nhap}  = 0;
            $ans->{$_}->{nq}   = $ans->{$_}->{nwithq} = $ans->{$_}->{qfreq} = 0;
            $ans->{$_}->{qci}  = [ (0, 0) ];
        }
        print "condense_sumstats:begin".Pretty(Dumper($ans)) if $o_verbose;
        my %region_seen;
        foreach (@{$sumstats->{$k}}) {
            my $rgn = pop_region($_->{pop});
            $region_seen{$rgn}++;
            my $nhap = $_->{nindiv} * $o_ploidy;
            my $nq   = $nhap - ($nhap * $_->{freqp});
            $ans->{$rgn}->{npop}   += 1;
            $ans->{$rgn}->{nindiv} += $_->{nindiv};
            $ans->{$rgn}->{nhap}   += $nhap;
            $ans->{$rgn}->{nq}     += $nq;
            $ans->{$rgn}->{nwithq} += ($_->{freqp} < 1.0 ? 1 : 0);
        }
        my $n_region_withq = 0;
        foreach (keys %region_seen) {  # region summaries and site totals
            $ans->{$_}->{qfreq} = $ans->{$_}->{nq} / $ans->{$_}->{nhap};
            $ans->{$_}->{qci}   = [ binomial_confint($ans->{$_}->{nhap}, $ans->{$_}->{nq}, $o_alpha, 0) ];
            # totals
            $ans->{_nindiv} += $ans->{$_}->{nindiv};
            $ans->{_nhap}   += $ans->{$_}->{nhap};
            $ans->{_nq}     += $ans->{$_}->{nq};
            $ans->{_nwithq} += $ans->{$_}->{nwithq};
            ++$n_region_withq if $ans->{$_}->{nwithq};
        }
        # totals
        $ans->{_qfreq} = $ans->{_nq} / $ans->{_nhap};
        $ans->{_qci}   = [ binomial_confint($ans->{_nhap}, $ans->{_nq}, $o_alpha, 0) ];
        $ans->{_region_seen} = \%region_seen;
        $ans->{_n_region} = scalar keys %region_seen;
        $ans->{_n_region_withq} = $n_region_withq;
        print "condense_sumstats:end:".Pretty(Dumper($ans)) if $o_verbose;
        $sumstats->{$k} = $ans;
    }
    # now rank the SNPs and add this rank to the locus
    # first criterion is highest _nwithq
    # second criterion is highest _n_region_withq
    # third criterion is highest _n_region
    # fourth criterion is highest lower-bound of _qci
    # fifth criterion is highest _qfreq
    # so the first is lowest rank, the last, highest
    my @k_rank = sort {
       $sumstats->{$a}->{_nwithq}         <=> $sumstats->{$b}->{_nwithq}         ||
       $sumstats->{$a}->{_n_region_withq} <=> $sumstats->{$b}->{_n_region_withq} ||
       $sumstats->{$a}->{_n_region}       <=> $sumstats->{$b}->{_n_region}       ||
       $sumstats->{$a}->{_qci}->[0]       <=> $sumstats->{$b}->{_qci}->[0]       ||
       $sumstats->{$a}->{_qfreq}          <=> $sumstats->{$b}->{_qfreq}
    } keys %$sumstats;
    foreach (0 .. $#k_rank) {
        $sumstats->{$k_rank[$_]}->{_rank} = $_ + 1;
    }
    print "condense_sumstats:veryend:".Pretty(Dumper($locus)) if $o_verbose;
}


sub binomial_confint($$$$) {
    # https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    # 4 digits right of decimal point
    # 1: number of trials
    # 2: number of successes
    # 3: alpha level
    # 4: method: 0: wilson  1: normal approx using adjusted p  2: normal approx
    my ($n, $success, $alpha, $method) = @_;
    my $z = Statistics::Distributions::udistr($alpha / 2);  # two-sized z-score
    usage("alpha $alpha seems inappropriate, z-score is $z") if $z < 0;
    my @ci;
    if (! $method) { # wilson's method
        my $z2 = $z * $z;
        my $plusminus = $z * sqrt(((1.0/$n) * $success * ($n - $success)) + 0.25*$z2);
        @ci = map { (1.0/($n + $z2)) * ($success + 0.5*$z2 + ($_ * $plusminus)) } (-1.0, +1.0);
    } elsif ($method == 1) { # adjust p and use normal approximation
        my $p = ($success + 2) / ($n + 4);
        my $plusminus = $z * sqrt((1.0 / $n) * $p * (1.0 - $p));  # the +/- part
        @ci = map { $p + ($_ * $plusminus) } (-1.0, +1.0);
    } elsif ($method == 2) { # normal approximation, alternative factorisation
        my $plusminus = $z * sqrt((1.0 / $n) * $success * ($n - $success));  # the +/- part
        @ci = map { (1.0 / $n) * ($success + ($_ * $plusminus)) } (-1.0, +1.0);
    } else { die "unknown method: $method"; }
    @ci = map { sprintf("%.${o_freqdigits}f", $_) } @ci;
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
        say Pretty(Dumper($locus)) if $o_verbose;
        my $lfs = $locus->{focal_sumstats};
        say STDOUT join("\t",
            $name,
            # _n_region(npop,npop,npop)
            "$lfs->{_n_region}(".
            "$lfs->{north}->{npop},".
            "$lfs->{other}->{npop},".
            "$lfs->{oland}->{npop})"
            ,
            # _n_region_withq(nwithq,nwithq,nwithq)
            "$lfs->{_n_region_withq}(".
            "$lfs->{north}->{nwithq},".
            "$lfs->{other}->{nwithq},".
            "$lfs->{oland}->{nwithq})"
            ,
            # _qfreq(qfreq,qfreq,qfreq)
            sprintf("%.${o_freqdigits}f", $lfs->{_qfreq})."(".
            sprintf("%.${o_freqdigits}f", $lfs->{north}->{qfreq}).",".
            sprintf("%.${o_freqdigits}f", $lfs->{other}->{qfreq}).",".
            sprintf("%.${o_freqdigits}f", $lfs->{oland}->{qfreq}).")"
            ,
            $locus->{focal_SNP} - 1,
            "$locus->{n_lsnps}:$locus->{n_rsnps}",
            $template);
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
    die "create_locus_template: locus $locus->{locus} has no focal SNP" if ! exists $locus->{focal_SNP};
    ## create template name
    my $name;
    if ($o_extended) {
        $name = $o_extended;
        # region coding
        my $R = $locus->{focal_sumstats}->{_n_region_withq};
        my ($n, $t, $l) = ($locus->{focal_sumstats}->{north}->{nwithq},
                           $locus->{focal_sumstats}->{other}->{nwithq},
                           $locus->{focal_sumstats}->{oland}->{nwithq});
        if ($R == 1) {
            $name .= "R1";
            $name .= ($n ? "N" : ($t ? "T" : "L"));
        } elsif ($R == 2) {
            $name .= "R2";
            $name .= ($n ? ($t ? "X" : "Y") : "Z");
        } elsif ($R == 3) {
            $name .= "R3A";
        } else {
            die "create_locus_template: locus $locus->{locus}: _n_region_withq makes no sense: $R";
        }
        # number of all populations with minor allele
        $name .= "P".sprintf("%.2d", $n + $t + $l);
        # populations with minor allele in each region
        $name .= "p$n$t$l";
        $name .= "_";
        # locus number
        $name .= sprintf("%.${o_padding}d", $locus->{locus});
    } else {
        $name = $o_name.sprintf("%.${o_padding}d", $locus->{locus});
    }
    ## create template sequence
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
        ++$N_locus_no_SNP;
        say STDOUT "no SNPs" if $o_verbose; 
        return 0;
    }

    # sort SNPs by column
    $locus->{snps} = [ sort { $a->{column} <=> $b->{column} } @{$locus->{snps}} ];
    ##my @s = sort { $a->{column} <=> $b->{column} } @{$LOCUS{$l}->{snps}};

    # reject templates with any 3+ allele SNPs
    my @snps3 = grep { $_->{al3} ne '-' } @{$locus->{snps}};
    if (@snps3) {
        ++$N_locus_SNP_3allele;
        say STDOUT 'at least one 3+ allele SNP' if $o_verbose; 
        return 0;
    }

    # reject template with no central SNP
    #my @csnps = grep { $_->{column} > $o_minflank and $_->{column} < ($locus->{seqlen} - $o_minflank + 1) } @{$locus->{snps}};
    my @csnps = grep { is_central_SNP($locus, $_->{column}) } @{$locus->{snps}};
    if (! @csnps) {
        ++$N_locus_no_central_SNP;
        say STDOUT 'no central snps' if $o_verbose;
        return 0;
    }
    my @i_csnps_has_sumstats = check_SNP_has_sumstats($locus, \@csnps);
    if (@i_csnps_has_sumstats == 0) {
        ++$N_locus_central_SNP_no_sumstats;
        say STDOUT "evaluate_locus: locus $locus->{locus}: no central SNP has sumstats" if $o_verbose;
        return 0;
    }

    # get SNPs on the flanks
    my @lsnps = grep { $_->{column} <= $o_minflank } @{$locus->{snps}};
    my $l_column = @lsnps > 0 ? $lsnps[$#lsnps]->{column} : 0;  # column of the rightmost left SNP
    my @rsnps = grep { $_->{column} >= ($locus->{seqlen} - $o_minflank + 1) } @{$locus->{snps}};
    my $r_column = @rsnps > 0 ? $rsnps[0]->{column} : 0;  # column of the leftmost right SNP

    my @i_csnps_ok = check_SNP_sumstats_OK($locus, \@csnps);
    if (@i_csnps_ok == 0) {
        ++$N_locus_central_SNP_no_sumstats_OK;
        say STDOUT "evaluate_locus: locus $locus->{locus}: no central SNP has sumstats OK" if $o_verbose;
        return 0;
    }

    if (@csnps == 1) {
        if (($l_column && ($csnps[0]->{column} - $l_column) < $o_centralsnpgap) or
            ($r_column && ($r_column - $csnps[0]->{column}) < $o_centralsnpgap)) {
            ++$N_locus_central_SNP_crowded;
            say STDOUT 'single central SNP too close to SNP in flank' if $o_verbose;
            return 0;
        }
    } else {  # more than one central SNP
        # none of the closeness checks require the SNPs to be sumstats_OK
        my @close = grep { ($csnps[$_ + 1]->{column} - $csnps[$_]->{column}) < $o_centralsnpgap }  0 .. ($#csnps - 1);
        if (@close) {
            ++$N_locus_central_SNP_crowded;
            say STDOUT 'at least two central SNPs are too close to each other' if $o_verbose;
            return 0;
        }
        if (($l_column && ($csnps[0]->{column} - $l_column) < $o_centralsnpgap) or
            ($r_column && ($r_column - $csnps[$#csnps]->{column}) < $o_centralsnpgap)) {
            ++$N_locus_central_SNP_crowded;
            say STDOUT 'at least one central SNP too close to SNP in flank' if $o_verbose;
            return 0;
        }
        if (@i_csnps_ok == 1) {  # if only one SNP OK, then focus on that
            ++$N_locus_central_SNP_multiple_one_OK;
            say STDOUT 'multiple central SNPs and only one is OK' if $o_verbose;
            while (@csnps > $i_csnps_ok[0] + 1) { # take non-OK SNPs off the right onto @rsnps
                unshift @rsnps, pop @csnps;
            }
            while (@csnps > 1) { # take non-OK SNPs off the left side onto @lsnps
                push @lsnps, shift @csnps;
            }
        } else {  # more than one OK SNP, pick the best by rank
            ++$N_locus_central_SNP_multiple_OK;
            # find the position in @csnps of the SNP with highest rank
            my ($i, $r) = (-1, 0); # i is the position, r is the highest rank seen so far
            foreach (@i_csnps_ok) {
                die "both i and r must be unset at same time" if $i == -1 xor $r == 0;
                if (($i == -1 && $r == 0) or $csnps[$_]->{sumstats}->{_rank} > $r) {
                    $i = $_;
                    $r = $csnps[$_]->{sumstats}->{_rank};
                }
            }
            say STDOUT "multiple central SNPs and the one at position $i with rank $r is the best" if $o_verbose;
            while (@csnps > $i + 1) { # take lower-rank SNPs off the right onto @rsnps
                unshift @rsnps, pop @csnps;
            }
            while (@csnps > 1) { # take lower-rank SNPs off the left side onto @lsnps
                push @lsnps, shift @csnps;
            }
        }
            
        # choose a central SNP, move the other(s) to the flanks
        #
        # 1) favour equal numbers of SNPs on both sides of central SNP
        # 2) favour a more central SNP for the focal SNP over a less central SNP
        while (@csnps > 1) {
            ++$N_locus_central_SNP_multiple_remaining;
            say STDOUT "*** After selecting among OK SNPs, it still seems \@csnps has more than 1 SNP ... this should not happen";
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
        ++$N_locus_flank_SNP_toomany;
        say STDOUT 'too many SNPs on left flank' if $o_verbose;
        return 0;
    }
    if (@rsnps > $o_maxflanksnps) {
        ++$N_locus_flank_SNP_toomany;
        say STDOUT 'too many SNPs on right flank' if $o_verbose;
        return 0;
    }

    ####### ####### ####### ####### ####### ####### #######

    say STDOUT "SUITABLE SNP at column $csnps[0]->{column}" if $o_verbose;

    $locus->{n_lsnps} = scalar @lsnps;
    $locus->{n_rsnps} = scalar @rsnps;
    $locus->{focal_SNP} = $csnps[0]->{column};
    $locus->{focal_sumstats} = $locus->{sumstats}->{$csnps[0]->{column}};
    return $locus->{focal_SNP};
}


sub check_SNP_has_sumstats($$) {
    my ($locus, $snps) = @_;
    my @ok;
    my $i = -1;
    foreach (@$snps) {
        ++$i;
        my $ok = (defined $locus->{sumstats}->{$_->{column}}) ? 1 : 0;
        push @ok, $i if $ok;
    }
    return @ok;
}


sub check_SNP_sumstats_OK($$) {
    my ($locus, $snps) = @_;
    my @ok;
    my $i = -1;
    foreach (@$snps) {
        ++$i;
        my $s = $locus->{sumstats}->{$_->{column}};
        my $ok = (defined $s and
               $s->{_n_region} >= $o_crit_n_region and
               $s->{_nwithq}   >= $o_crit_nwithq   and
               $s->{_qfreq}    >= $o_crit_qfreq)       ? 1 : 0;
        $_->{sumstats_OK} = $ok;
        $_->{sumstats} = $s;
        push @ok, $i if $ok;
    }
    return @ok;
}


