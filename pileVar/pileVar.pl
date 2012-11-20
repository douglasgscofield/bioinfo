#! /usr/bin/perl -w

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
#   read_bases_column = 4
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
use POSIX;
use Getopt::Long;
use Pod::Usage;

my $in_file;
my $seq_column = 0;
my $coord_column = 1;
my $ref_base_column = 2;
my $cvrg_column = 3;
my $read_bases_column = 4;
my $base_quality_column = 5;
my $quality_cutoff = 0; # phred scale integer
my $mincov = 0;
my $minqual = 0;
my $minqualcov = 0;
my $variants_only = 0;
my $bed = "No"; #set to "Yes" to convert coordinates to bed format (0-based start, 1-based end)
my $out_file; #  = $ARGV[10];
my $print_bases_quals = 0;
my $no_coord_check = 0;

my $stdio;
my $help;
my $noheader;
my $ploidy = 1;
my $consensus = 1;
my $consensus_seq_name = "consensus_seq";
my $consensus_seq = "";
my $consensus_fasta = "";
my $fasta_gap_char = "n";
my $qual = 1;
my $qualcalc = "rms";
my $qualround = 1;
my $base_qual_offset = 64; # by default, assume Phred+64 quality scores
my $track_indels = 0;
my $fixed_frac = 0.8;
my $filter_indel_frac = 0.05;

my $usage = "
NAME

  pileVar.pl - parse SAMtools mpileup format for variant information

SYNOPSIS

  pileVar.pl  [options]

OPTIONS

    -                        read input from STDIN, write to STDOUT
    -in <filename>           read input from <filename>, else from STDIN
    -out <filename>          write output to <filename>, else to STDOUT
    -noheader                do not print header on output
    -base_qual_offset <int>  offset of base quality ASCII value from 0 quality [default $base_qual_offset]

    -ploidy <int>            ploidy of samples from which pileup is derived [default $ploidy]

    -no_coord_check          do not check coordinates for increase-by-1 consistency.  If a consensus
                             FASTA file is being produced, this suppresses the insertion of gaps
                             in the output
    -variants_only           print out only those coordinates containing variants from the reference
                             [default $variants_only]
    -consensus               call consensus sequence from reads, adds two columns to output,
                             one for consensus call, one for mean quality of call [default $consensus]
    -consensus_fasta <filename>  print the consensus sequence in FASTA format to the given file
    -fasta_gap_char <char>   the character used for coordinate-skipping gaps in the consensus
                             FASTA file, gaps are also identified by position in the FASTA header line 
                             [default $fasta_gap_char]
    -print_bases_quals       print out bases and base qualities from input mpileup
                             [default $print_bases_quals]

  Quality:

    -mincov <int>            minimum raw coverage to call a change to reference [default $mincov]
    -minqual <int>           minimum base quality cutoff [default $minqual]
    -minqualcov <int>        minimum coverage of quality bases to call a change to reference [default $minqualcov]
    -qual                    print qualities of variants [default $qual]
    -qualcalc mean|rms|sum   method for calculating qualities of variants [default $qualcalc]
    -qualround <int>         digit to which to round variant qualities [default $qualround]

  Indels:

    -track_indels            track the presence of indels [default $track_indels]
    -filter_indel_frac <fraction>  do not report indels at a coordinate if the fraction of 
                             reads containing them is below <fraction> [default $filter_indel_frac]

  SNP variants:

    -fixed_frac <fraction>   minimum fraction of non-ref bases required to call a base as 
                             fixed with respect to the reference [default $fixed_frac]

    -help, -?                help message

";

sub print_usage_and_exit($) {
    my $msg = shift;
    print "$msg\n" if $msg;
    print $usage;
    exit 0;
}

GetOptions(
    "" => \$stdio, 
    "in=s" => \$in_file,
    "out=s" => \$out_file,
    "noheader" => \$noheader,
    "base_qual_offset" => \$base_qual_offset,
    "ploidy=i" => \$ploidy,
    "no_coord_check" => \$no_coord_check,
    "variants_only:i" => \$variants_only,
    "consensus:i" => \$consensus,
    "consensus_fasta=s" => \$consensus_fasta,
    "fasta_gap_char=s" => \$fasta_gap_char,
    "print_bases_quals" => \$print_bases_quals,
    "mincov=i" => \$mincov,
    "minqual=i" => \$minqual,
    "minqualcov=i" => \$minqualcov,
    "qual" => \$qual,
    "qualcalc=s" => \$qualcalc,
    "qualround=i" => \$qualround,
    "track_indels" => \$track_indels,
    "filter_indel_frac=f" => \$filter_indel_frac,
    "fixed_frac=f" => \$fixed_frac,
    "help|?" => \$help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $help;
print_usage_and_exit("invalid value for option -qualcalc '$qualcalc'") if $qualcalc ne "mean" 
    and $qualcalc ne "rms"
    and $qualcalc ne "sum";
print_usage_and_exit("invalid value for option -qualround '$qualround', must be >= 0") if $qualround < 0;
print_usage_and_exit("invalid value for option -fixed_frac '$fixed_frac', must be between 0 and 1") if $fixed_frac < 0 or $fixed_frac > 1;
print_usage_and_exit("invalid value for option -filter_indel_frac '$filter_indel_frac' must be between 0 and 1") if $filter_indel_frac < 0 or $filter_indel_frac > 1;

my $invalid_line_counter = 0;
my $first_skipped_line = "";

if ($stdio or ! defined($in_file)) {
    *IN = *STDIN;
} else {
    open (IN, "<$in_file") or die "Cannot open $in_file $!\n";
}

if ($stdio or ! defined($out_file)) {
    *OUT = *STDOUT;
} else {
    open (OUT, ">$out_file") or die "Cannot open $out_file $!\n";
}

if (! $noheader) {
    print OUT "sequence\t";
    print OUT "coord\t";
    print OUT "ref\t";
    print OUT "cvg\t";
    print OUT "bases\tquals\t" if $print_bases_quals;
    print OUT "A\tC\tG\tT\t";
    print OUT "qualA\tqualC\tqualG\tqualT\t" if $qual;
    print OUT "qual_cvg\t";
    print OUT "ref_diff\t";
    print OUT "cons\tcons_qual\t" if $consensus;
    print OUT "note\n";
}

my $prev_coord = -1;

while (<IN>) {
    chomp;
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
    my $indel_in_reads = 0;
    my $indel_del_reads = 0;
    my $mincov_not_met = 0;
    my $minqualcov_not_met = 0;
    my $note = "";
    if ($fields[ $ref_base_column ] eq "*") {
        print "skipping indicated indel line '*'\n" if $track_indels;
        next; # skip indel lines
    }

    if (! $no_coord_check and $prev_coord >= 0 and ($prev_coord + 1) != $fields[ $coord_column ]) {
        my $coord_gap = $fields[$coord_column] - $prev_coord - 1;
	$note .= ", " if $note;
	$note .= "coord $fields[$coord_column] out of sequence with prev $prev_coord, gap of $coord_gap";
        if ($consensus_fasta and $coord_gap > 0) {  # insert a gap, if it is a gap
            $consensus_seq .= $fasta_gap_char x $coord_gap;
            $consensus_seq_name .= "|gap[" . ($prev_coord+1) . "-" . ($fields[$coord_column]-1) . "]";
        }
	print STDERR "coord $fields[$coord_column] out of sequence with prev $prev_coord, gap of $coord_gap\n";
    }

    my $read_bases   = $fields[ $read_bases_column ];
    die "Coverage column" . ($cvrg_column+1) . " contains non-numeric values. Check your input parameters as well as format of input dataset." if ( not isdigit $fields[ $cvrg_column ] );
    $mincov_not_met = 1 if $fields[ $cvrg_column ] < $mincov;
    my $base_quality = $fields[ $base_quality_column ];
    if ($read_bases =~ m/[\$\^\+-]/) {
        $read_bases =~ s/\^.//g; #removing the start of the read segement mark
        $read_bases =~ s/\$//g; #removing end of the read segment mark
        while ($read_bases =~ m/([\+-]){1}(\d+)/g) {
            my $indel_type = $1;
            my $indel_len = $2;
            $read_bases =~ s/([\+-]{1}$indel_len.{$indel_len})//; # remove indel info from read base field
            my $indel_contents = $1;
            ++$indels{"indel, $indel_type, $indel_len bases, $indel_contents"};
            ++$indel_in_reads if $indel_type eq "+";
            ++$indel_del_reads if $indel_type eq "-";
        }
    }

    my $N_bases = length($read_bases);

    my $indel_frac_coord = ($indel_in_reads + $indel_del_reads) / $N_bases;

    if ( length($read_bases) != length($base_quality) ) {
        $first_skipped_line = $. if $first_skipped_line eq "";
        ++$invalid_line_counter;
        next;
    }
    # after removing read block and indel data the length of read_base 
    # field should identical to the length of base_quality field
    
    my @bases = split //, $read_bases;
    my @qv    = split //, $base_quality;

    for my $base ( 0 .. @bases - 1 ) {
        my $this_base_qual = ord( $qv[ $base ] ) - $base_qual_offset;
        if ( $this_base_qual >= $minqual and $bases[ $base ] ne '*')
        {
            ++$above_qv_bases;
            
            if ( $bases[ $base ] =~ m/[ATGC]/i )
            {
                $SNPs_exist = 1;    
                $SNPs{ uc( $bases[ $base ] ) } += 1;
                $SNPqualsum{ uc( $bases[ $base ] ) } += $this_base_qual;
                $SNPqualsumsq{ uc( $bases[ $base ] ) } += ($this_base_qual * $this_base_qual);
                $diff_count += 1;
            } elsif ( $bases[ $base ] =~ m/[\.,]/ ) {
                $SNPs{ uc( $fields[ $ref_base_column ] ) } += 1;
                $SNPqualsum{ uc( $fields[ $ref_base_column ] ) } += $this_base_qual;
                $SNPqualsumsq{ uc( $fields[ $ref_base_column ] ) } += ($this_base_qual * $this_base_qual);
            }           
        }
    } 
    $minqualcov_not_met = 1 if $above_qv_bases < $minqualcov;
    
    if ($bed eq "Yes") {
           my $start = $fields[ $coord_column ] - 1;
           my $end   = $fields[ $coord_column ];
           $fields[ $coord_column ] = "$start\t$end";
    } 
    
    push @output_fields, $fields[ $seq_column ];
    push @output_fields, $fields[ $coord_column ];
    push @output_fields, $fields[ $ref_base_column ];
    push @output_fields, $fields[ $cvrg_column ];
    push @output_fields, $fields[ $read_bases_column ] if $print_bases_quals;
    push @output_fields, $fields[ $base_quality_column ] if $print_bases_quals;

    foreach my $SNP (sort keys %SNPs) {
        push @output_fields, $SNPs{$SNP};
    }

    foreach my $SNP (sort keys %SNPs) {
        my $qualval;

        if ($SNPs{$SNP} == 0) {
           $qualval = 0;
        } elsif ($qualcalc eq "sum") {
           $qualval = $SNPqualsum{$SNP};
        } elsif ($qualcalc eq "mean") {
           $qualval = $SNPqualsum{$SNP} / $SNPs{$SNP};
        } elsif ($qualcalc eq "rms") {
           $qualval = sqrt($SNPqualsumsq{$SNP} / $SNPs{$SNP});
        }
        $SNPqual{$SNP} = $qualval;
        $qualval = sprintf "%.${qualround}f", $qualval if $qualval != 0;

        push @output_fields, $qualval if $qual;
    }
    
    push @output_fields, $above_qv_bases;
    push @output_fields, $diff_count;

    my $ref_base = uc( $fields[ $ref_base_column ] );
    my $consensus_base = $ref_base;
    my $consensus_qual =  $SNPqual{$ref_base};

    if ($mincov_not_met or $minqualcov_not_met) {

        # do not call a change to the reference

        if ($mincov_not_met) {
	    $note .= ", " if $note;
            $note .= "below min coverage $mincov";
        }
        if ($minqualcov_not_met) {
	    $note .= ", " if $note;
            $note .= "below min qual coverage $minqualcov";
        }

    } else {

        # check to see if consensus different from reference

        foreach my $SNP (sort keys %SNPs) {
                my $SNP_frac = $above_qv_bases ? ($SNPs{$SNP} / $above_qv_bases) : 0;
                if ($SNP_frac >= $fixed_frac and $SNP ne $ref_base) {
                    $consensus_base = $SNP;
                    $consensus_qual = $SNPqual{$SNP};
	            $note .= ", " if $note;
                    $note .= ("fixed $ref_base -> $SNP (" . sprintf("%.3f", $SNP_frac) . ")");
                }
        }
    }

    # notes for each coordinate

    if ($indel_frac_coord >= $filter_indel_frac) {
	$note .= ", " if $note;
        $note .= "indel fraction";
        $note .= " +" . sprintf("%.3f", $indel_in_reads / $N_bases) if ($indel_in_reads / $N_bases) > 0;
        $note .= " -" . sprintf("%.3f", $indel_del_reads / $N_bases) if ($indel_del_reads / $N_bases) > 0;
    }

    # consensus columns

    if ($consensus) {
	push @output_fields, $consensus_base;
	push @output_fields, sprintf("%.${qualround}f", $consensus_qual);
    }
    $consensus_seq .= $consensus_base if $consensus_fasta;

    push @output_fields, $note if $note;

    print OUT join("\t", @output_fields), "\n";

    if ($track_indels and $indel_frac_coord >= $filter_indel_frac) {
        foreach my $indel_msg (sort keys %indels) {
            print $indel_msg, ", $indels{$indel_msg} times\n"; 
        }
    }

    # end of iteration stuff

    $prev_coord = $fields[$coord_column];
    print STDERR "coord $prev_coord\r" if ! ($prev_coord % 100);
}

print "Skipped $invalid_line_counter invalid line(s) beginning with line $first_skipped_line\n" if $invalid_line_counter > 0;

if ($consensus_fasta) {
    # I could redo this to peel off FASTA lines as sufficient length of consensus
    # sequence is produced.  It could save a lot of memory.
    open (FASTA, ">$consensus_fasta") or die "Cannot open $consensus_fasta $!\n";
    print FASTA ">", $consensus_seq_name, "|length=", length($consensus_seq), "\n";
    my $fasta_line_len = 72;
    while (my $line = substr($consensus_seq, 0, $fasta_line_len, "")) {
        print FASTA "$line\n";
    }
    close FASTA;
}

