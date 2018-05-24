#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use Data::Dumper;

my $debug = 0;

use constant missing_rho => -1;
use constant missing_rho_output => 'NA';

my $o_rhocolumn = 4;
my $o_rhofile = "";
my $o_faifile = "";
my $o_posfile = "";
my $o_bedfile = "";
my %FAI_LEN;  # lengths of known reference sequences, loaded from $o_faifile
my $o_help;

sub usage {
    say STDERR join("\n", @_) if $#_;
    say STDERR "

$0 - calculate rho for intervals between SNP positions

Infer recombination rates (rho) in arbitrary genomic intervals, using
length-weighted contributions within underlying intervals.  Three files are
expected for input, all with tab-separated columns:

* a four-column BED file describing recombination rate estimates in genomic intervals (--rhofile)
* a file specifying positions between which recombination rates should be inferred.  This file 
  can have one of two formats, determined by the option used to specify the file
  - a two-column reference,position file describing genomic positions between which values of rho
    should be inferred (--posfile).  If this file is provided, rho is inferred for the intervals
    between positions, as well as the first and last intervals of each chromosome, if they are
    not covered by positions here.  **OR**
  - a three-column BED file describing genomic intervals within which rho should be inferred (--bedfile).
    If this file is provided, rho is inferred for only these intervals.
* an FAI file describing the sizes of reference sequences named within the other files (--faifile)


Options
-------

    --rhofile FILE    BED file of rho intervals
    --faifile FILE    File containing Fasta index (.fai) for reference

    --posfile FILE    File with SNP ref and position described in cols 1 and 2  **OR**
    --bedfile FILE    BED file specifying genomic intervals

    --rho-column INT  Column of --rhofile containing the rho estimate, numbered from 1 [default 4]

    --help | -?       help

";
    exit 1;
}

GetOptions( "rhofile=s"    => \$o_rhofile,
            "faifile=s"    => \$o_faifile,
            "posfile=s"    => \$o_posfile,
            "bedfile=s"    => \$o_bedfile,
            "rho-column=i" => \$o_rhocolumn,
            "help|?"       => \$o_help
) or usage();

usage() if $o_help or ! $o_rhofile or ! $o_faifile or ! ($o_posfile xor $o_bedfile);

sub load_rho_intervals($$);
sub create_rho_ref_hash($);
sub calc_rho($$$$);
sub rho_interval_from_ref_pos($);
sub rho_interval_from_bed($);

my ($RI, $RH) = load_rho_intervals($o_rhofile, $o_faifile);  # now also returns reference to hash

if ($o_posfile) {
    rho_interval_from_ref_pos($o_posfile);
} elsif ($o_bedfile) {
    rho_interval_from_bed($o_bedfile);
} else {
    usage("oops, didn't provide --posfile or --bedfile");
}


# this is called if --bedfile was specified
sub rho_interval_from_ref_pos($) {
    # read SNPs and dummy up a BED interval (but +1 left side) with ref prev-snp this-snp rho
    my $file = shift;
    my ($prev_ref, $prev_pos) = ("", 0);
    open (my $f, "<", $file) or die "could not open --posfile '$file': $!";
    while (<$f>) {
        chomp;
        my ($ref, $pos) = split/\t/;
        $pos -= 1;  # make it 0-based
        die "reference name not found in FAI $ref" if not defined $FAI_LEN{$ref};
        if ($prev_ref eq "" or $prev_ref ne $ref) {
            if ($prev_ref ne "") { # finish off the previous chromosome
                my $this_rho = calc_rho($prev_ref, $prev_pos, $prev_ref, $FAI_LEN{$prev_ref});
                say STDOUT join("\t", $prev_ref, $prev_pos, $FAI_LEN{$ref}, $this_rho);
            }
            ($prev_ref, $prev_pos) = ($ref, 0);
        }
        die "$file:$.: consecutive positions out of order" if $pos <= $prev_pos;
        my $this_rho = calc_rho($prev_ref, $prev_pos, $ref, $pos);
        say STDOUT join("\t", $ref, $prev_pos, $pos, $this_rho);
        ($prev_ref, $prev_pos) = ($ref, $pos);
    }
}

# this is called if --posfile was specified
sub rho_interval_from_bed($) {
    # read BED intervals and calculate rho for each
    my $file = shift;
    open (my $f, "<", $file) or die "could not open --bedfile '$file': $!";
    while (<$f>) {
        chomp;
        my ($ref, $left, $right) = split/\t/;
        die "reference name not found in FAI $ref" if not defined $FAI_LEN{$ref};
        my $this_rho = calc_rho($ref, $left, $ref, $right);
        say STDOUT join("\t", $ref, $left, $right, $this_rho);
    }
}

#
# Subroutines to load rho data, create rho interval data structures, and
# calculate arbitrary rho intervals
#
#----------------------------------

sub calc_rho($$$$) {
    my ($ref1, $pos1, $ref2, $pos2) = @_;
    my $debug = 0;
    if ((! $ref1 and ! $pos1) and (! $ref2 and ! $pos2)) {
        die "calc_rho: both ref/pos pairs not set";
    } elsif ((! $ref1 and ! $pos1) xor (! $ref2 and ! $pos2)) {
        say STDERR "$0:calc_rho: one ref/pos pair not set, returning ".missing_rho_output;
        return missing_rho_output;
    }
    my $rho;
    if ($ref1 eq $ref2) {
        die "calc_rho: position 1 $pos1 >= position 2 $pos2" if $pos1 >= $pos2;
        my $i1; # interval for $ref1/$pos1
        my $i2; # interval for $ref2/$pos2
        my $is_NA = 0; # count of intervals with missing rho values
        for ($i1 = $RH->{$ref1};    ! ($pos1 >= $RI->[$i1]->[1] and $pos1 < $RI->[$i1]->[2]); ++$i1 ) { }
        say STDERR "calc_rho: 1: $ref1:$pos1 found within interval $i1: ".join(" ", @{$RI->[$i1]}) if $debug > 1;
        # now step through intervals until we find the complete list spanning snp1 and snp2, keeping track
        # if any intervening intervals are missing rho estimates
        for ($i2 = $i1; ! ($pos2 > $RI->[$i2]->[1] and $pos2 <= $RI->[$i2]->[2]); ++$i2 ) {
            say STDERR "rho for $i2 is $RI->[$i2]->[3]" if $debug > 1;
            ++$is_NA if $RI->[$i2]->[3] == missing_rho;
        }
        ++$is_NA if $RI->[$i2]->[3] == missing_rho;
        say STDERR "calc_rho: 2: $ref2:$pos2 found within interval $i2: ".join(" ", @{$RI->[$i2]}) if $debug > 1;
        say STDERR "calc_rho: $is_NA intervals were missing rho values" if $debug > 1;

        return missing_rho_output if $is_NA;  # stop trying to calculate rho if there are missing values

        my $rho_sum = 0;
        foreach ($i1 .. $i2) {
            my $i_bp = $RI->[$_]->[2] - $RI->[$_]->[1];
            my $left_crop  = $pos1 - $RI->[$_]->[1]; # bases to substract b/c of pos1
            my $right_crop = $RI->[$_]->[2] - $pos2; # bases to substract b/c of pos2
            $i_bp -= $left_crop  if $left_crop > 0;
            $i_bp -= $right_crop if $right_crop > 0;
            $rho_sum += ($i_bp * $RI->[$_]->[3]);
            say STDERR "interval $_: [$RI->[$_]->[1] $RI->[$_]->[2] rho=$RI->[$_]->[3]] i_bp=$i_bp pos1=$pos1 pos2=$pos2 left_crop=$left_crop right_crop=$right_crop final i_bp=$i_bp rho_sum=$rho_sum" if $debug > 1;
        }
        $rho = $rho_sum / ($pos2 - $pos1);
        say STDERR "final rho from $pos1 to $pos2: $rho" if $debug > 1;

    } else {
        return missing_rho_output;  # stop if positions are on different refs
    }
    return $rho;
}

#----------------------------------

sub load_rho_intervals($$) {
    my ($o_rhofile, $o_faifile) = @_;
    my $prev_chr = "";
    my $prev_pos = 0;
    my $debug = 0;
    my @RI;
    say STDERR "$0:load_rho_intervals: assuming proper fai file is $o_faifile" if $debug;
    open (my $fai, "<", $o_faifile) or die "could not open fai file $o_faifile $!";
    while (<$fai>) {
        my @l = split/\t/;
        $FAI_LEN{$l[0]} = $l[1];
    }
    close($fai);
    say STDERR "$0:load_rho_intervals: read ".scalar(keys %FAI_LEN)." fai records" if $debug;

    say STDERR "$0:load_rho_intervals: attempting to open recombination map file '$o_rhofile'" if $debug;

    open (my $rhofd, "<", $o_rhofile) or die "could not open rho file $o_rhofile $!";
    while (<$rhofd>) {
        chomp;
        my @l = split/\t/;
        die "$o_rhofile$.: interval endpoints out of order" if $l[2] <= $l[1];
        my $this_rho = missing_rho;
        $this_rho = $l[$o_rhocolumn - 1] if $l[$o_rhocolumn - 1] ne 'NA';
        die "$o_rhofile$.: rho value '$this_rho' 0 or unknown" if $this_rho == 0;
        die "$o_rhofile$.: chromosome name makes no sense" if $l[0] eq "";
        if ($l[0] ne $prev_chr) { # new or first chromosome
            if ($prev_chr ne "") {  # if new chromosome (NOT first), did we miss any at the end?
                my $fai_end = $FAI_LEN{$prev_chr};
                die "unknown reference name in FAI $prev_chr" if not defined $fai_end;
                if ($prev_pos < $fai_end) {
                    say STDERR "$o_rhofile$.:load_rho_intervals: $l[0]: no rho at end of $prev_chr from pos $prev_pos to $fai_end" if $debug > 1;
                    # create a gap to cover the missing region at the end of the chromosome
                    push @RI, [$prev_chr, $prev_pos, $fai_end, missing_rho];
                }
            }
            $prev_chr = $l[0];
            $prev_pos = 0;
        }
        if ($prev_pos != $l[1]) {
            say STDERR "$o_rhofile$.:load_rho_intervals: $l[0]: no rho from pos $prev_pos to $l[1]" if $debug > 1;
            # create a gap to cover the missing region before processing the current gap
            push @RI, [$l[0], $prev_pos, $l[1], missing_rho];
        }
        # process the current gap
        push @RI, [$l[0], $l[1], $l[2], $this_rho];
        $prev_pos = $l[2];
    }
    my $fai_end = $FAI_LEN{$prev_chr};
    die "unknown reference name in FAI $prev_chr" if not defined $fai_end;
    if ($prev_pos < $fai_end) {
        say STDERR "$o_rhofile$.:load_rho_intervals: $prev_chr: end of $o_rhofile no rho at end of $prev_chr from pos $prev_pos to $fai_end" if $debug > 1;
        # create a gap to cover the missing region at the end of the chromosome
        push @RI, [$prev_chr, $prev_pos, $fai_end, missing_rho];
    }

    if ($debug > 1) {
        foreach (0 .. $#RI) {
            say STDERR "RI $_: ".join('  ', @{$RI[$_]});
            if ($_ > 0 && $RI[$_]->[0] eq $RI[$_ - 1]->[0] and $RI[$_]->[1] != $RI[$_ - 1]->[2]) {
                say STDERR "**** BREAK in intervals at $_";
            }
        }
    }

    my $RH = create_rho_ref_hash(\@RI);

    return (\@RI, $RH);
}

#----------------------------------

# hash with keys holding indices to first interval for reference in @{$RI}

sub create_rho_ref_hash($) {
    my $RI = shift;
    my %RH;
    my $debug = 0;
    my $prev_chr = "";
    say STDERR "I think there are 0 .. ".(scalar(@{$RI})-1)." entries in \$RI" if $debug > 1;
    foreach (0..(scalar(@{$RI})-1)) {
        if ($RI->[$_]->[0] ne $prev_chr) {
            die "create_rho_ref_hash: adding $RI->[$_]->[0] but it is already defined" if defined $RH{$RI->[$_]->[0]};
            $RH{$RI->[$_]->[0]} = $_;
            $prev_chr = $RI->[$_]->[0];
        }
    }
    return \%RH;
}

