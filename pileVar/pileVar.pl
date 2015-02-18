#!/usr/bin/env perl

# Call variants based on output of samtools mpileup command.
#
# Built upon the structure of the very useful Galaxy mpileup-format parser script
# https://bitbucket.org/galaxy/galaxy-central/src/tip/tools/samtools/pileup_parser.pl.
#
# Douglas G. Scofield, Umeå Plant Sciences Centre, Umeå Sweden
# douglas.scofield@plantphys.umu.se.  
# douglasgscofield@gmail.com
#
# Input is **raw** mpileup output, via something like
#
#      samtools mpileup -BQ0 -d 1000000 -f $Reference -D $BAM -A
#
# This requests no read filtering, and includes anomalous base pairs (-A).  You might want to
# drop that last one, based on the settings of BAM flags following your pipeline.
#
# The initial development of this script is for the calling of variants during 
# reference-guided assembly of a chloroplast sequence.  Thus, my first interest
# is uncovering *fixed differences* between the reference and the set of reads,
# rather than finding high-confidence SNPs/SNVs that are easy to isolate against
# a background of sequence variation.  Eventually I will extend this or write 
# another script to handle SNP calling in pooled samples.
#
#
# CHANGELOG
#
# 2012-11-29
# - option name changes to bring them to my standard practices
# - change naming of --no-pos-check option
# - implement --variants-only option
# 2012-02-27
# - change flag variable $cvrg_cutoff to $mincov to match the command-line option
# - change flag variable $quality_cutoff to $minqual to match the command-line option
# - add flag -minqualcov and flag $minqualcov to impose a cutoff on the number of bases 
#   that match the quality criterion
# - for coordinates that do not meet $mincov or $minqualcov, create a full output line with 
#   the base kept as the reference, and a note indicating that minimum coverage was not met.
#   Previously the script would simply skip such positions
#
# 2012-02-23
# - get rid of the original pileup_parser.pl's use of $out_string, use @output_fields
#
# 2012-02-22
# - print gap in output FASTA if gap detected in coordinates, not done with -no_coord_check
#
# 2012-02-21
# - print warning if coordinates are not monotonically increasing by 1 each line, add option
#   -no_coord_check to suppress this
#
# 2012-02-20
# - add option -filter_indel_frac <fraction> which does not report indels at a coordinate when
#   the fraction of reads involved is below the given fraction
# - add -consensus option to produce column of consensus sequence, currently only pays attention
#   to -fixed_frac option
# - change $SNPs_only to $variants_only
# - various streamlinings of code
#
# 2012-02-17
# - assign defaults for many of the columns based on samtools' mpileup format (as of samtools 0.1.18)
#   coord_column = 1
#   ref_base_column = 2
#   cvrg_column = 3
#   bas_col = 4
#   base_quality_column = 5
# - assign defaults for script behaviors
#   quality_cutoff = 0, do not assign a quality cutoff when counting bases
#   cvrg_cutoff = 0, do not exclude bases below a particular quality
#   SNPs_only = "No", prints everything, including non-SNP coordinates
#   bed = "No", do not convert coordinates to BED format
#   total_diff removed, always print the column of total differences from the reference
#   print_qual_bases = "No", do not print the base and quality columns
# - use Getopt::Long and provide the following options
#   -                      read from STDIN, write to STDOUT
#   -in <string>           read from file, else STDIN
#   -out <string>          write to file, else STDOUT
#   -noheader              do not print header on output
#   -ploidy <int>          assume pileup is derived from a sample having the given ploidy
#   -fixed_frac <fraction> the minimum fraction of quality bases differing from the 
#                          reference that are required to call a position as fixed
#   -mincov <int>          minimum coverage cutoff 
#   -minqual <int>         minimum quality cutoff to include a base
#   -qual                  provide quality columns for each base
#   -qualcalc <string>     method for calculating quality for each variant
#   -qualround <int>       digit to which to round quality value
#   -base_qual_offset <int>  offset of base quality
#   -track_indels          track presence of indels
#   -help|-?               help message
# - indel recognition with -track_indels, print a message about each indel and its
#   appearance in reads
# - fixation of changes with -fixed_frac, print a note about each SNP above a given
#   fractional representation, ex: "fixed A -> G (0.994)"
#
# TODO
# problem calling heterozygotes!
# MA_1    4667    C       85      .,t,,..TT.,T.t,.tT$tt,TttT.TT,.T.TT.TT.,,tTT,.,,tT,t.T.TTt,.,TTtt,T,.,.,,,,TT.tTt.TT,,  ce.ifg[46fi,h0ae.!..c0.05f..fa.e0.e!!efe04.gef^5.h:e!f/,/f^c,.1.e.eihiihfe04i,..g=-ff
# MA_1    5528    T       77      .CCc,cC,C.C.c.,c^].C.,..,,C.,C,.CCCC,C..cCC.cCC,CC..C,.ccc,.C.,C,,.,C,Ccc.,.,,C R!3/g.4i1f.f4fd<E2ddeedf1fe2dc1111`;Rd/22c.21i41df.[d3.4ef2cf2iihc1h/0/hfadi/
# 
#
# - handle mincov violations not by skipping, but by noting min coverage was not met
#   and keeping reference base
# - continue improving SNP quality calculations
# - recognize beginnings/ends of reads and base memberships in reads
# - recognize read mapping qualities
# - handle excluding bases belonging to reads mapping below a certain quality
# - allow an option to rescale base quality based on both mapping and call quality,
#   like MAQ and others but with my own flavor
# - improve indel reporting by tracking "net" indel depth across coordinates
# 

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use List::Util qw/sum/;

my $in_file;
my $seq_col = 0;
my $pos_col = 1;
my $ref_col = 2;
my $cvg_col = 3;
my $bas_col = 4;
my $qua_col = 5;
my $map_col = 6;
my $opt_mappingquality = 0;
my $opt_mincov = 0;
my $opt_minqual = 0;
my $opt_minqualcov = 0;
my $opt_variantsonly = 0;
my $opt_hetz = 0;
my $opt_hetzminfreq = 0.1;
my $opt_hetzcheck = 1;
my $out_file; #  = $ARGV[10];
my $opt_printbasesquals = 0;
my $opt_noposcheck = 0;

my $opt_stdio;
my $opt_help;
my $opt_noheader;
my $opt_ignoreN = 0;
my $opt_ploidy = 2;
my $opt_consensus = 0;
my $consensus_seq_name = "consensus_seq";
my $consensus_seq = "";
my $opt_consensusfasta = "";
my $opt_fastagapchar = "n";
my $opt_qual = 0;
my $opt_qualcalc = "rms";
my $opt_qualround = 1;
my $opt_basequaloffset = 64; # by default, assume Phred+64 quality scores
my $opt_indels = 0;
my $opt_indelmode = 0;
my $opt_fixedfrac = 0.8;
my $opt_indelfrac = 0.05;

my $usage = "
NAME

  pileVar.pl - parse SAMtools mpileup format for variant information

SYNOPSIS

  pileVar.pl  [options]

OPTIONS

    -                        read input from stdin, write to stdout
    -i FILE, --input FILE    read input from FILE, else from stdin
    -o FILE, --output FILE   write output to FILE, else to stdout
    --mapping-quality        input has mapping quality column (samtools mpileup -s)
    --no-header              do not print header on output
    --base-qual-offset INT   offset of base quality ASCII value from 0 quality [default $opt_basequaloffset]

    --ignore-N               do not produce output for positions with reference of N

    --ploidy INT             ploidy of sample, currently ignored [default $opt_ploidy]

    --no-pos-check           do not check positions for increase-by-1 consistency.  If a
                             consensus FASTA file is being produced, this suppresses the 
                             insertion of gaps in the output (--no-coord-check is a synonym)
    --hetz                   determine heterozygous positions, applying all other options
    --hetz-min-freq FLOAT    minimum frequency for heterozygous call at a site [default $opt_hetzminfreq]
    --hetz-check             check that #alleles do not exceed ploidy, if so ignore it [default $opt_hetzcheck]
                             this all would be better replaced with a model
    --variants-only          print out only those positions containing variants from the 
                             reference [default $opt_variantsonly]
    --consensus              call consensus sequence from reads, adds two columns to output,
                             one for consensus call, one for mean quality of call [default $opt_consensus]
    --consensus-fasta FILE   print the consensus sequence in FASTA format to the given file
    --fasta-gap-char CHAR    the character used for position-skipping gaps in the consensus
                             FASTA file, gaps are also identified by position in the FASTA 
                             header line [default $opt_fastagapchar]
    --print-bases-quals      print out bases and base qualities from input mpileup
                             [default $opt_printbasesquals]

  Quality:

    --mincov INT             minimum raw coverage to call a change to reference [default $opt_mincov]
    --minqual INT            minimum base quality cutoff [default $opt_minqual]
    --minqualcov INT         minimum coverage of quality bases to call a change to 
                             reference [default $opt_minqualcov]
    --qual                   print qualities of variants [default $opt_qual]
    --qualcalc mean|rms|sum  method for calculating qualities of variants [default $opt_qualcalc]
    --qualround INT          digit to which to round variant qualities [default $opt_qualround]

  Indels:

    --indels                 track the presence of indels [default $opt_indels]
    --indel-frac FLOAT       do not report indels at a position if the fraction of 
                             reads containing them is below FLOAT [default $opt_indelfrac]
    --indel-mode             track ONLY the presence of indels [default $opt_indelmode]

  SNP variants:

    --fixed-frac FLOAT       minimum fraction of non-ref bases required to call a base as 
                             fixed with respect to the reference [default $opt_fixedfrac]

    --help, -?               help message

";

sub print_usage_and_exit($) {
    my $msg = shift;
    print "$msg\n" if $msg;
    print $usage;
    exit 0;
}

GetOptions(
    "" => \$opt_stdio, 
    "input|i=s" => \$in_file,
    "output|o=s" => \$out_file,
    "mapping-quality" => \$opt_mappingquality,
    "no-header" => \$opt_noheader,
    "base-qual-offset" => \$opt_basequaloffset,
    "ignore-N" => \$opt_ignoreN,
    "ploidy=i" => \$opt_ploidy,
    "no-pos-check|no-coord-check" => \$opt_noposcheck,
    "hetz" => \$opt_hetz,
    "hetz-min-freq" => \$opt_hetzminfreq,
    "hetz-check" => \$opt_hetzcheck,
    "variants-only" => \$opt_variantsonly,
    "consensus" => \$opt_consensus,
    "consensus-fasta=s" => \$opt_consensusfasta,
    "fasta-gap-char=s" => \$opt_fastagapchar,
    "print-bases-quals" => \$opt_printbasesquals,
    "mincov=i" => \$opt_mincov,
    "minqual=i" => \$opt_minqual,
    "minqualcov=i" => \$opt_minqualcov,
    "qual" => \$opt_qual,
    "qualcalc=s" => \$opt_qualcalc,
    "qualround=i" => \$opt_qualround,
    "indels" => \$opt_indels,
    "indel-frac=f" => \$opt_indelfrac,
    "indel-mode" => \$opt_indelmode,
    "fixed-frac=f" => \$opt_fixedfrac,
    "help|?" => \$opt_help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $opt_help;
print_usage_and_exit("invalid value for option --qualcalc '$opt_qualcalc'") if $opt_qualcalc ne "mean" 
    and $opt_qualcalc ne "rms"
    and $opt_qualcalc ne "sum";
print_usage_and_exit("invalid value for option --qualround '$opt_qualround', must be >= 0") if $opt_qualround < 0;
print_usage_and_exit("invalid value for option --fixed-frac '$opt_fixedfrac', must be between 0 and 1") if $opt_fixedfrac < 0 or $opt_fixedfrac > 1;
print_usage_and_exit("invalid value for option --indel-frac '$opt_indelfrac' must be between 0 and 1") if $opt_indelfrac < 0 or $opt_indelfrac > 1;
print_usage_and_exit("--hetz must be specified with --ploidy of >= 2") if $opt_hetz and $opt_ploidy == 1;
print_usage_and_exit("invalid value for option --hetz-min-freq '$opt_hetzminfreq' must be greater than 0 and less than 1") if $opt_hetzminfreq <= 0 or $opt_hetzminfreq >= 1;

my $invalid_line_counter = 0;
my $first_skipped_line = "";

if ($opt_stdio or ! defined($in_file)) {
    *IN = *STDIN;
} else {
    open (IN, "<$in_file") or die "Cannot open $in_file $!\n";
}

if ($opt_stdio or ! defined($out_file)) {
    *OUT = *STDOUT;
} else {
    open (OUT, ">$out_file") or die "Cannot open $out_file $!\n";
}

if (! $opt_noheader) {
    print OUT "seq\t";
    print OUT "pos\t";
    print OUT "ref\t";
    print OUT "cvg\t";
    if ($opt_indelmode) {
        print OUT "iqual\tifreq\tilen\tioper\tn_var\n";
    } else {
        print OUT "bases\tquals\t" if $opt_printbasesquals;
        if ($opt_hetz) {
            print OUT "gt\tfreq\tqual\t";
        } else {
	    print OUT "A\tC\tG\tT\t";
        }
        print OUT "qualA\tqualC\tqualG\tqualT\t" if $opt_qual;
        print OUT "minqual\t" if $opt_minqual;
        print OUT "n_var\t";
        print OUT "cons\tcons_qual\t" if $opt_consensus;
        print OUT "note\n";
    }
}

my $prev_pos = -1;

while (<IN>) {
    chomp;
    my $entire_line = $_;
    next if m/^\#/;
    my @fields = split /\t/;
    my @output_fields = ();
    use constant SNPs_zeroes => ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0);
    my %SNPs = SNPs_zeroes;
    my %SNPqual = SNPs_zeroes;
    my %SNPqualsum = SNPs_zeroes;
    my %SNPqualsumsq = SNPs_zeroes;
    my $above_qv_bases = 0;
    my $SNPs_exist = 0;
    my $diff_count = 0;
    my %indels = ();
    my %indel_operations = ();
    my %indel_lengths = ();
    my $indel_in_reads = 0;
    my $indel_del_reads = 0;
    my $mincov_not_met = 0;
    my $minqualcov_not_met = 0;
    my $note = "";
    if ($fields[ $ref_col ] eq "*") {
        print "skipping deletion in reference '*'\n"; # if $opt_indels;
        next; # skip indel lines
    }
    if ($opt_ignoreN and uc($fields[ $ref_col ]) eq "N") {
        next;
    }

    my $reference_sequence_name = $fields[ $seq_col ];
    push @output_fields, $reference_sequence_name;
    push @output_fields, $fields[ $pos_col ];
    push @output_fields, $fields[ $ref_col ];
    push @output_fields, $fields[ $cvg_col ];
    push @output_fields, $fields[ $bas_col ] if $opt_printbasesquals;
    push @output_fields, $fields[ $qua_col ] if $opt_printbasesquals;

    if (! $opt_noposcheck and $prev_pos >= 0 and ($prev_pos + 1) != $fields[ $pos_col ]) {
        my $pos_gap = $fields[$pos_col] - $prev_pos - 1;
	$note .= ", " if $note;
	$note .= "pos $fields[$pos_col] out of sequence with prev $prev_pos, gap of $pos_gap";
        if ($opt_consensusfasta and $pos_gap > 0) {  # insert a gap, if it is a gap
            $consensus_seq .= $opt_fastagapchar x $pos_gap;
            $consensus_seq_name .= "|gap[" . ($prev_pos+1) . "-" . ($fields[$pos_col]-1) . "]";
        }
	print STDERR "pos $fields[$pos_col] out of sequence with prev $prev_pos, gap of $pos_gap\n";
    }

    my $read_bases   = $fields[ $bas_col ];
    die "Coverage column" . ($cvg_col+1) . " contains non-numeric values. Check your input parameters as well as format of input dataset." if ( not isdigit $fields[ $cvg_col ] );
    $mincov_not_met = 1 if $fields[ $cvg_col ] < $opt_mincov;
    my $base_quality = $fields[ $qua_col ];
    if ($read_bases =~ m/[\$\^\+-]/) {
        $read_bases =~ s/\^.//g; #removing the start of the read segement mark
        $read_bases =~ s/\$//g; #removing end of the read segment mark
        while ($read_bases =~ m/([\+-]){1}(\d+)/g) {
            my $indel_type = $1;
            my $indel_len = $2;
            if ($indel_len > 1024) { # there is a problem with this indel field and probably this whole pileup line...
                print STDERR "line $. skipped, invalid indel len $indel_len: line is '$entire_line'\n";
                $first_skipped_line = $. if $first_skipped_line eq "";
                ++$invalid_line_counter;
                next;
            }
            $read_bases =~ s/([\+-]{1}$indel_len.{$indel_len})//; # remove indel info from read base field
            my $indel_contents = $1;
            ++$indels{"indel, $indel_type, $indel_len bases, $indel_contents"};
            ++$indel_lengths{($indel_type eq "+" ? +$indel_len : -$indel_len)};
            ++$indel_operations{uc($indel_contents)};
            ++$indel_in_reads if $indel_type eq "+";
            ++$indel_del_reads if $indel_type eq "-";
        }
    }

    my $N_bases = length($read_bases);

    my $indel_frac_pos = ($indel_in_reads + $indel_del_reads) / $N_bases;

    if ($opt_indelmode) {
        next if ! $indel_frac_pos;  # non-indels not needed
        # @output_fields already has $reference_sequence_name,
        # $fields[ $pos_col ], $fields[ $ref_col ], $fields[ $cvg_col ]
        if ($indel_frac_pos >= $opt_indelfrac and scalar(keys %indel_lengths) == 1 
            and scalar(keys %indel_operations) == 1) {
            # above cutoff frequency and just one length and operation
            push @output_fields, "good";
        } else {
            push @output_fields, "bad";
        }
        push @output_fields, (sprintf "%.4f", $indel_frac_pos);
        push @output_fields, join(",", keys %indel_lengths);
        push @output_fields, join(",", keys %indel_operations);
        push @output_fields, $indel_in_reads + $indel_del_reads;
    	print OUT join("\t", @output_fields), "\n";
        next;  # no more is needed
    }

    if ( length($read_bases) != length($base_quality) ) {
        $first_skipped_line = $. if $first_skipped_line eq "";
        ++$invalid_line_counter;
        next;
    }
    # after removing read block and indel data the length of read_base 
    # field should identical to the length of base_quality field
    
    my @bases = split //, $read_bases;
    my @qv    = split //, $base_quality;

    for my $i ( 0 .. @bases - 1 ) {

        my $this_base_qual = ord( $qv[ $i ] ) - $opt_basequaloffset;

        if ( ! $opt_minqual or ($this_base_qual >= $opt_minqual and $bases[ $i ] ne '*')) {

            ++$above_qv_bases;
            
            if ( $bases[ $i ] =~ m/[ATGC]/i ) {
                $SNPs_exist = 1;    
                $SNPs{ uc( $bases[ $i ] ) } += 1;
                $SNPqualsum{ uc( $bases[ $i ] ) } += $this_base_qual;
                $SNPqualsumsq{ uc( $bases[ $i ] ) } += ($this_base_qual * $this_base_qual);
                $diff_count += 1;
            } elsif ( $bases[ $i ] =~ m/[\.,]/ ) {
                $SNPs{ uc( $fields[ $ref_col ] ) } += 1;
                $SNPqualsum{ uc( $fields[ $ref_col ] ) } += $this_base_qual;
                $SNPqualsumsq{ uc( $fields[ $ref_col ] ) } += ($this_base_qual * $this_base_qual);
            }           

        }
    } 

    $minqualcov_not_met = 1 if $above_qv_bases < $opt_minqualcov;
    
    if ($opt_hetz) { # print gt, minor-allele freq, qual columns

        my $minor_allele = "";
        my $minor_allele_freq = 1.1;  # always set the first time through
        my @gt = ();  # genotypes other than the minor allele
        my $gt_qual;  # simply look at the qualities of the bases, not ideal by any means

        foreach my $SNP (sort keys %SNPs) {
            my $SNP_freq = $above_qv_bases ? ($SNPs{$SNP} / $above_qv_bases) : 0;
            if ($SNP_freq >= $opt_hetzminfreq) {  # it's acceptable
                if ($SNP_freq < $minor_allele_freq) { # switch places
                    push @gt, $minor_allele if $minor_allele ne "";
                    $minor_allele = $SNP;
                    $minor_allele_freq = $SNP_freq;
                } else {
                    push @gt, $SNP;
                }
            }
        }
        push @gt, $minor_allele;  # add minor allele to end

        if ($above_qv_bases == 0 or $minqualcov_not_met or $mincov_not_met or ($opt_hetzcheck and scalar(@gt) > $opt_ploidy)) {
            @gt = ( $fields[ $ref_col ] );  # force to reference, not ideal by any means
            # still show frequency to be true frequency of reference allele
            $minor_allele_freq = $above_qv_bases ? ($SNPs{$minor_allele} / $above_qv_bases) : 0;
            $gt_qual = $above_qv_bases ? ($SNPs{$minor_allele} / $above_qv_bases) : 0;
        } else {
            if ($opt_qualcalc eq "sum") {
               $gt_qual = sum(@SNPqualsum{@gt});
            } elsif ($opt_qualcalc eq "mean") {
               $gt_qual = sum(@SNPqualsum{@gt}) / sum(@SNPs{@gt});
            } elsif ($opt_qualcalc eq "rms") {
               $gt_qual = sqrt(sum(@SNPqualsumsq{@gt}) / sum(@SNPs{@gt}));
            }
        }

        push @output_fields, join("/", @gt);
        push @output_fields, (sprintf "%0.4f", $minor_allele_freq);
        push @output_fields, (sprintf "%.${opt_qualround}f", $gt_qual);

    } else { # print base counts instead

        foreach my $SNP (sort keys %SNPs) {
            push @output_fields, $SNPs{$SNP};
        }

    }

    foreach my $SNP (sort keys %SNPs) {
        my $qualval;

        if ($SNPs{$SNP} == 0) {
           $qualval = 0;
        } elsif ($opt_qualcalc eq "sum") {
           $qualval = $SNPqualsum{$SNP};
        } elsif ($opt_qualcalc eq "mean") {
           $qualval = $SNPqualsum{$SNP} / $SNPs{$SNP};
        } elsif ($opt_qualcalc eq "rms") {
           $qualval = sqrt($SNPqualsumsq{$SNP} / $SNPs{$SNP});
        }
        $SNPqual{$SNP} = $qualval;
        $qualval = sprintf "%.${opt_qualround}f", $qualval if $qualval != 0;

        push @output_fields, $qualval if $opt_qual;
    }
    
    push @output_fields, $above_qv_bases if $opt_minqual;
    push @output_fields, $diff_count;

    my $ref_base = uc( $fields[ $ref_col ] );
    my $consensus_base = $ref_base;
    my $consensus_qual =  $SNPqual{$ref_base};

    if ($mincov_not_met or $minqualcov_not_met) {

        # do not call a change to the reference

        if ($mincov_not_met) {
	    $note .= ", " if $note;
            $note .= "below min coverage $opt_mincov";
        }
        if ($minqualcov_not_met) {
	    $note .= ", " if $note;
            $note .= "below min qual coverage $opt_minqualcov";
        }

    } else {

        # check to see if consensus different from reference

        foreach my $SNP (sort keys %SNPs) {
            my $SNP_frac = $above_qv_bases ? ($SNPs{$SNP} / $above_qv_bases) : 0;
            if ($SNP_frac >= $opt_fixedfrac and $SNP ne $ref_base) {
                $consensus_base = $SNP;
                $consensus_qual = $SNPqual{$SNP};
	        $note .= ", " if $note;
                $note .= ("fixed $ref_base -> $SNP (" . sprintf("%.3f", $SNP_frac) . ")");
            }
        }
    }

    # notes for each position

    if ($indel_frac_pos >= $opt_indelfrac) {
	$note .= ", " if $note;
        $note .= "indel fraction";
        $note .= " +" . sprintf("%.3f", $indel_in_reads / $N_bases) if ($indel_in_reads / $N_bases) > 0;
        $note .= " -" . sprintf("%.3f", $indel_del_reads / $N_bases) if ($indel_del_reads / $N_bases) > 0;
    }

    # consensus columns

    if ($opt_consensus) {
	push @output_fields, $consensus_base;
	push @output_fields, sprintf("%.${opt_qualround}f", $consensus_qual);
    }
    $consensus_seq .= $consensus_base if $opt_consensusfasta;

    push @output_fields, $note if $note;

    if (! $opt_variantsonly || ($opt_variantsonly && ($SNPs_exist || $indel_in_reads || $indel_del_reads))) {
    	print OUT join("\t", @output_fields), "\n";

    	if ($opt_indels and $indel_frac_pos >= $opt_indelfrac) {
            foreach my $indel_msg (sort keys %indels) {
                print $indel_msg, ", $indels{$indel_msg} times\n"; 
	    }
    	}
    }

    # end of iteration stuff

    $prev_pos = $fields[$pos_col];
    print STDERR "ref $reference_sequence_name pos $prev_pos        \r" if ! ($prev_pos % 100);
}

print "Skipped $invalid_line_counter invalid line(s) beginning with line $first_skipped_line\n" if $invalid_line_counter > 0;

if ($opt_consensusfasta) {
    # I could redo this to peel off FASTA lines as sufficient length of consensus
    # sequence is produced.  It could save a lot of memory.
    open (FASTA, ">$opt_consensusfasta") or die "Cannot open $opt_consensusfasta $!\n";
    print FASTA ">", $consensus_seq_name, "|length=", length($consensus_seq), "\n";
    my $fasta_line_len = 72;
    while (my $line = substr($consensus_seq, 0, $fasta_line_len, "")) {
        print FASTA "$line\n";
    }
    close FASTA;
}

