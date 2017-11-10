Bioinformatics scripts
======================

These have been useful to me, I hope they can be useful to you!


plink-pairwise-loci.pl
----------------------

    USAGE:   plink-pairwise-loci.pl [ --help ] --ped file.ped --mibs file.mibs

Scan `file.ped` and write matrix to STDOUT counting the number of complete
pairwise comparisons, with no missing alleles, between each genotype described
in `file.ped`.  The `file.mibs` will be used to verify the format of the output but
will not otherwise be used.

A locus with missing alleles is indicated by `0 0` in `file.ped`.  A message
is produced if only one allele is coded as missing, and the entire locus is
coded as missing.

`file.ped` is formatted for PLINK input and `file.mibs` is PLINK output of
pairwise inbreeding coefficients, which is produced without regard to the
number of complete pairwise comparisons.


fai2WindowBed
-------------

Perl script to create evenly-spaced BED intervals along Fasta sequences described by the given Fasta index file. A Fasta index file can be produced using '`samtools faidx ...`' and other tools. This is useful for creating BED files to guide sliding-window analyses with other tools, for example '`bedtools intersect`'.  If the `--overlap` option is non-zero, windows will overlap by the given number of bp.

    SYNOPSIS

        scripts/fai2WindowBed.pl --fai file.fa.fai --width 100000 > file.windows100000.bed

    Create BED file of evenly-spaced windows from Fasta index.

    This script creates evenly-spaced windows along sequences described in a Fasta
    index file ('.fai').  It is useful for creating BED files for sliding window
    analyses in other tools.

    Window width is set with --width, and window overlap (if any) is set with
    --overlap.  The BED intervals are named after the sequence and the coordinate
    of the beginning of the window plus 1 to match the convention in non-BED files,
    i.e.,

          chr01 \t 0 \t 100000 \t chr01_00000001
          chr01 \t 100000 \t 200000 \t chr01_00100001

    The last window for each sequence is shortened, if necessary, not to extend
    beyond the end of the sequence; underwidth windows can be dropped from the
    output by specifying the --drop-underwidth option.

    The BED file of windows is written to stdout with sequences appearing in their
    order within the '.fai' file and BED intervals in increasing order of starting
    coordinate.

    OPTIONS

        --fai FILE         REQUIRED: Fasta index file; this can be created with
                           'samtools faidx file.fa' and other tools
        --width INT        Window width [10000]
        --overlap INT      Window overlap [0]
        --no-pad           Do not add 0-padding to the digits used to write the
                           starting coordinate in the interval name; without this,
                           the number of digits is calculated from the maximum
                           needed to make all coordinates equal-width in the
                           intervals across each sequence.
        --no-add-one       Do not add 1 to the interval name
        --drop-underwidth  Do not include underwidth windows in the output; this
                           drops the last window for each sequence if it is not
                           the specified width.



fillBed
-------

Perl script to fill BED intervals of Fasta sequences with a given single-character sequence.  Especially useful for converting BED intervals to `N`.  The case of the substituted sequence is conserved by default, so an interval containing `TTatccgC` would be substituted (using `-s N`) as `NNnnnnnN`.  Requires BioPerl.

    SYNOPSIS

        fillBed -i input.fa -b to_be_filled.bed -s N -p -o output.fa

    This script does no validation of the sequence provided with -s, but BioPerl
    will do simple validation of the alphabet when updating the modified sequence.

    OPTIONS

        -i FILE      input FASTA sequences (read from STDIN if not specified)
        -o FILE      output FASTA sequences (written to STDOUT if not specified)
        -b FILE      BED-format file containing intervals
        -s CHAR      single character to fill intervals
        -p           preserve case of replaced characters [default 1]
        -P           DO NOT preserve case of replaced characters
        -?, --help   help message



fastaSort, gffSort
------------------

`fastaSort`: Sort Fasta file by ID, naturally (id_1, id_2, ..., id_10, id_11, ...).  One argument, the Fasta filename, and writes sorted output to standard output.  Requires `BioPerl 1.6.922` or thereabouts (`Bio::DB::Fasta` needs to have the `get_all_primary_ids` method), as well as `Sort::Naturally`.  Because of its use of `Bio::DB::Fasta`, it creates a specialised index file for the Fasta file, which can be time-consuming to create but greatly speeds up sorting once completed.  If you would like to remove this index after completion of the script, set `$remove_index` to 1 within the script.

`gffSort`: Sort GFF file by sequence name (column 1) then numerically by position (column 4).  One argument, the GFF filename, and writes sorted output to standard output while also removing comment lines specified with `###`, which are inserted by MAKER.


stacksExtractStats.pl
---------------------

Extract sample-specific stack and coverage information from a Stacks log file.  For the following input file `ex.log`, part of a Stacks log file:

    Identifying unique stacks; file   1 of  40 [XYZ_H0122_GATTAC.1.cat]
    /xxxx/yyyy/private/GBS_working_data/stacks-1.20/ustacks -t fastq -f ./XYZ/concat/XYZ_H0122_GATTAC.1.cat.fq -o ./XYZ_read1 -i 1 -m 3 -M 3 -p 15 -d -r 2>&1
    Min depth of coverage to create a stack: 3
    Max distance allowed between stacks: 3
    Max distance allowed to align secondary reads: 5
    Max number of stacks allowed per de novo locus: 3
    Deleveraging algorithm: enabled
    Removal algorithm: enabled
    Model type: SNP
    Alpha significance level for model: 0.05
    Parsing ./XYZ/concat/XYZ_H0122_GATTAC.1.cat.fq
    Loaded 2800 RAD-Tags; inserted 2423 elements into the RAD-Tags hash map.
      0 reads contained uncalled nucleotides that were modified.
      Mean coverage depth is 9; Std Dev: 18.0875 Max: 87
    Coverage mean: 9; stdev: 18.0875
    Deleveraging trigger: 27; Removal trigger: 45
    Calculating distance for removing repetitive stacks.
      Distance allowed between stacks: 1
      Using a k-mer length of 71
      Number of kmers per sequence: 74
      Miniumum number of k-mers to define a match: 3
    Removing repetitive stacks.
      Removed 3 stacks.
      38 stacks remain for merging.
    Calculating distance between stacks...
      Distance allowed between stacks: 3
      Using a k-mer length of 35
      Number of kmers per sequence: 110
      Miniumum number of k-mers to define a match: 5
    Merging stacks, maximum allowed distance: 3 nucleotide(s)
      38 stacks merged into 36 stacks; deleveraged 0 stacks; removed 0 stacks.
      Mean merged coverage depth is 9.41667; Std Dev: 20.8037; Max: 96
    Merging remainder radtags
      2461 remainder sequences left to merge.
      Distance allowed between stacks: 5
      Using a k-mer length of 23
      Number of kmers per sequence: 122
      Miniumum number of k-mers to define a match: 7
      Matched 29 remainder reads; unable to match 2432 remainder reads.
    Number of utilized reads: 368
    Writing loci, SNPs, and alleles to './XYZ_read1/'...
      Refetching sequencing IDs from ./XYZ/concat/XYZ_H0122_GATTAC.1.cat.fq... read 2800 sequence IDs.
    done.
   
The command

    ./stacksExtractStats.pl ex.log

produces a tab-separated table containing a header line and columns for sample name (the pattern recognized can easily be modified), sample number in the run, total samples in the run, number of RAD-Tags, number of stacks identified, and mean merged coverage depth.

    samplename        nsample totsamples      nradtags        nstacks coverage
    XYZ_H0122_GATTAC.1        1       40      2800    36      9.41667


mummer2Vcf.pl
-------------

Convert a list of SNPs and indels produced by [Mummer's](http://mummer.sourceforge.net/) `show-snps -T` command to a rough substitute for a [VCF file](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41).  VCF requires the first reference base to be output for indels and Mummer doesn't output this, so this script requires a file containing the Fasta-format reference sequence, in SNP order.  It will open the file named as the reference in the first line of the `show-snps -T`-formatted file, but this can also be supplied with the `--fasta` option, which will override the file-encoded reference sequence if it is provided.

It also provides the `--type SNP|INDEL` option to subselect the set of variants produced if you like.

Note the order of the SNPs and indels with respect to their reference sequences must be in the same order as the sequences within the Fasta file supplied with `--fasta`.  Otherwise, `mummer2Vcf.pl` will produce a '`couldn't find reference`' error.

The 'VCF' format produced by this script is approaching actual VCF, thanks to some very helpful user contributions.

This script requires BioPerl.

subsampleReads.pl
-----------------

    subsampleReads.pl  [ options ] [ - | fastq-file fastq-file ]

Input files are assumed to be paired-end, interleaved FastQ. Single-end input
data is specified with --single.  Input files may be gzipped (with suffix
'.gz').

    -c|--count INT        number of reads to keep... OR
    -f|--fraction FLOAT   fraction of reads to keep; a read pair is
                          selected if a random uniform draw <= FLOAT
    -l|--limit INT        with --fraction, take no more than INT reads, otherwise
                          ignored
    -s|--single           reads are single-end
    -                     read input from STDIN, write to STDOUT


To handle paired-end reads in separate files for reads 1 and 2, use this
script in a pipe like

```bash
shuffleFastq.pl - read1.fq.gz read2.fq.gz \
| subsampleReads.pl - -f 0.01 \
| deshuffleFastq.pl - sub1.fq.gz sub2.fq.gz
```


seqStats.pl
-----------

    seqStats.pl [ --interleaved ] [ - | fastq_file.fq[.gz] ]

Produce very simple output of FastQ-format sequence statistics: the number of sequences, the total base count also in Mbp, and the mean read length.  If the `--interleaved` option is provided, it also prints out an `npairs` line that is calculated simply by dividing the number of FastQ sequences by 2.

```bash
$ seqStats.pl file.1.fq.gz
nseq	975215
nbase	98496715
nMbp	98.496715
meanlen	101
```

phredDetector.pl
----------------

Heuristically determine the [Phred-scaled quality score](http://en.wikipedia.org/wiki/FASTQ_format) used in the FastQ file presented on input.  Can handle gzipped input.  If everything looks reasonable, it simply prints `33`, `64` or `59` (this last for older Solexa sequences with the `--solexa` option) to `stdout`.  It can be used in pipeline scripts to autodetect Phred encodings:

```bash
Quality=$(phredDetector.pl --solexa fastq_file.fq.gz)
if [ "$Quality" = "33" ] ; then
   echo "perhaps Illumina 1.8+"
elif [ "$Quality" = "64" ] ; then
   echo "perhaps Illumina 1.3+ to pre-1.8"
elif [ "$Quality" = "59" ] ; then
   echo "perhaps Solexa"
else
   echo "Couldn't autodetect quality, return value was '$Quality'"
fi
```


BioPerl is not required.  Usage is:

    phredDetector.pl  [ options ]  fastq_file1.fq[.gz]

Input must be in FastQ format, if a filename is given it may be gzipped (*.gz)

Output to stdout is either

* `33`, appears to be Illumina 1.8+ or Sanger quality encoding
* `64`, appears to be Illumina 1.3+ to pre-1.8 quality coding
* `59`, appears to be Solexa (pre-Illumina) base-64 quality (only with `--solexa`)
* `??`, the input couldn't be interpreted

The only data interpreted are read quality scores on line 4 of each read; all sequence, pairing information etc. is ignored.

**Caveat:** If all bases on input are of unusually high quality, then a Phred base of 64 (or 59 with `--solexa`)
may be reported when a Phred scale based on 33 was the one actually used.  A
few heuristics are used to detect possible problems, but these are not
comprehensive.

* If (maximum quality &minus; minimum quality) &ge; 47, a warning message is printed
  to stderr and detected quality encoding (possibly erroneous) to stdout; this value can be adjusted with the `--wide` option
* If (maximum quality &minus; minimum quality) &le; 10, a warning message is printed
  to stderr and detected quality encoding (possibly erroneous) to stdout; this value can be adjusted with the `--narrow` option.  If you run this script on quality-trimmed reads, you may trigger this warning.
* If otherwise unusual quality scores or unknown input were detected, an error
  message is printed to stderr and '??' to stdout

    --solexa
This script does not diagnose faulty FastQ files, nor does it fully diagnose
the various versions of Solexa, Sanger, Illumina pipelines described in
http://en.wikipedia.org/wiki/FASTQ_format.  It simply applies a decreasing minimum cutoff
from 64, to 59 (with `--solexa`), to 33.  It is thus compatible with Sanger encoding but
it does not take into account finer distinctions such as Illumina 1.3+ pipelines
not typically producing quality scores 0 and 1.

**Options**

    -             Read uncompressed FastQ from stdin
    --reads INT   Number of reads to process to determine Phred basis [10000]
                  If 0, process *all* reads in the input file
    --wide INT    Use INT for the 'too wide' first heuristic above [47]
    --narrow INT  Use INT for the 'too narrow' second heuristic above [10]
    --solexa      Add Solexa base-64 (beginning with 59 ';') as one of the types
                  to guess.  If this type of quality encoding is detected, '59'
                  is produced on stdout.  Without this option, only '33', '64' or 
                  error '??' are reported.
    --no-solexa   Do not include Solexa base-64 ('59') as one of the quality types
                  to guess.  This is the default.
    --verbose     Output includes the input filename followed by a TAB character

    --help | -?   Generate this help output


fermiExtractContigs.pl
------------------------

Extract Fasta-format contigs from Fermi's FastQ-like output files.  Fermi is an
overlap assembler developed by Heng Li (https://github.com/lh3/fermi).  

Fermi produces two formats of FastQ-like output: scaftig sequences in files
named *.fq.gz, and unitig sequences in files named *.mag.gz which use an
overlap graph file format called MAG.  This script can read both formats, with
scaftigs being the default input and MAG format interpreted with the `--mag` or
`--mag-verbose` options.

Output of this script is annotated, line-wrapped Fasta format, to stdout.  For
scaftig files, the FastQ-like read name line includes two additional
tab-separated fields encoding the length of the scaftig and the number of
non-redundant reads that built the scaftig.  This script includes these in the
Fasta description, and also calculates the median non-redundant read coverage
along the contig by interpreting the encoded coverage in fourth line.  So this:

````
@26417937:25351227_0	191	83
TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTATATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTAAACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCCTGCATGTAGAA
+
"#$%%&'((()**+,-../0123456789:;<=>?@ABCDEFGHIJKLLMNOPQRSTUVWXYZ[\]^_`abcdeffghijklmnoopqrstttttttttttsrqpponmmmlkkjihggfedcba`_^]\[ZYXWVUTSRQPONMLKJIIHGFEDCBA@?>=<;:9876543210//.-,+*)('&&%$#"
````

becomes this:

```
>26417937:25351227_0 length:191;n_reads:83;median_coverage:44
TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTAT
ATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTA
AACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCC
TGCATGTAGAA
```

[MAG-format input](https://github.com/lh3/fermi/blob/master/fermi.1) presents different fields in the read name line.  In MAG, the
second field is the number of non-redundant reads, and the following two fields
list the left and right neighbors of the unitig in the overlap graph structure.
Unitigs may have multiple left or right neighbours, or none, in which case the
neighbour information is a single '.'.  Because unitig length is not encoded in
the header, it is calculated from each unitig.  In the example below, 13 reads
were used to build the unitig, and it has two neighbour unitigs to the left,
one with 77 and one with 76 bases of exact overlap, and four neighbor unitigs
to the right, with 67, 60, 68 and 69 bases of exact overlap.

````
@9012595342:2229517731	13	355005708,77;958384802,76;	2817374678,67;3657728097,60;7691722363,68;9232081066,69; 
TGTTTTATTAATAAAATATCTCTCTAACTTGTTTATTTCTGACCTGTTTCAGGTGAATTCGAAATTGCATGGGATTTGATGGATAGTTTAGAGGATATTTGGATGATTTATTTTCATTTTAGTTTCCTAGTTTAGCTGATCTTGGGAATATCTTCCCAGATTATGTAAACAGTTTGATAACTTCCACAGGGAGGTTTTACCCTGTGGAAAACTTATAAATACTTATTAT
+
"""""""""##########$$$$$%%&&&&&&&&&&&&&&&&&&&&&&&&&'''''''''(((((((((()))))))*******+++++++++++++++++*********)))))))***)))))))())))))))))))))))))))))))(((((((((''''''''''&&&&&&&%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$########"""
````

Output with `--mag` is

```
>9012595342:2229517731 length:229;n_reads:13;median_coverage:6
TGTTTTATTAATAAAATATCTCTCTAACTTGTTTATTTCTGACCTGTTTCAGGTGAATTC
GAAATTGCATGGGATTTGATGGATAGTTTAGAGGATATTTGGATGATTTATTTTCATTTT
AGTTTCCTAGTTTAGCTGATCTTGGGAATATCTTCCCAGATTATGTAAACAGTTTGATAA
CTTCCACAGGGAGGTTTTACCCTGTGGAAAACTTATAAATACTTATTAT
```

Output with `--mag-verbose` includes the lists of left and right neighbours:

```
>9012595342:2229517731 length:229;n_reads:13;median_coverage:6 355005708,77;958384802,76; 2817374678,67;3657728097,60;7691722363,68;9232081066,69;
TGTTTTATTAATAAAATATCTCTCTAACTTGTTTATTTCTGACCTGTTTCAGGTGAATTC
GAAATTGCATGGGATTTGATGGATAGTTTAGAGGATATTTGGATGATTTATTTTCATTTT
AGTTTCCTAGTTTAGCTGATCTTGGGAATATCTTCCCAGATTATGTAAACAGTTTGATAA
CTTCCACAGGGAGGTTTTACCCTGTGGAAAACTTATAAATACTTATTAT
```

This script required BioPerl to compose and write the Fasta-format output.



fermiExtractContigs_simple.sh
------------------------

Extract Fasta-format contigs from a
[fermi](https://github.com/lh3/fermi)-format `*.fq.gz` FastQ-like scaftig
files.  Does not use Perl.  Writes Fasta to `stdout`, with each sequence given its fermi sequence name and a description that includes sequence length and number of non-redundant reads that built the scaftig.  So this:

```
@26417937:25351227_0	191	83
TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTATATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTAAACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCCTGCATGTAGAA
+
"#$%%&'((()**+,-../0123456789:;<=>?@ABCDEFGHIJKLLMNOPQRSTUVWXYZ[\]^_`abcdeffghijklmnoopqrstttttttttttsrqpponmmmlkkjihggfedcba`_^]\[ZYXWVUTSRQPONMLKJIIHGFEDCBA@?>=<;:9876543210//.-,+*)('&&%$#"
```

Becomes this:

```
>26417937:25351227_0 length:191,n_reads:83
TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTATATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTAAACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCCTGCATGTAGAA
```

**Differences from `fermiExtractContigs.pl`**: In the interests of speed and simplicity, this script does no wrapping of the fermi-generated sequences, likewise it does not include median coverage along the entire length of the sequence, nor can it extract the proper fields from fermi's MAG-format headers.  If these limitations are not suitable for your task, consider using the `fermiExtractContigs.pl` script above.


windowWig
---------

Say you have a data stream (optionally containing a header line) with the 
following format:

    reference1  1        value1
    reference1  2        value2
    ...
    reference1  1100102  value3
    reference2  1        value4
    reference2  2        value5
    ...

Such data might represent, say, mapped-read coverage values by position
within reference sequences.  This script will create a USCS [WIG][] `fixedStep`
file by summarizing values within fixed-size windows across positions within
references.  It currently summarizes values by their median (see below for
more details) but this can easily be changed.  Also, there are some 
common-sense ordering requirements on the input data, see below.

As a quick example, say you want to produce a WIG file summarizing median
coverage across a reference genome, you want to use all read mappings, and
the mappings are contained within a few separate BAM files.  A fast pipeline
for doing this could be:

```bash
samtools mpileup -AB -d1000000 -q0 -Q0 -f ref.fa *.bam \
| mergePileupColumns \
| cut -f1,2,4 \
| windowWig > cov.wig
```

See below for `mergePileupColumns`.  Output in the WIG file will look something like:

    fixedStep chrom=reference1 start=1 step=50 span=50
    9
    19
    17
    15
    16
    23
    ...

This script walks the data stream, crossing each reference in column 1
in windowsize (default 50) chunks based on the position in the 2nd column, 
computing the median of the values in column 3 within each window.  The 
median is calculated based on the number of values seen within each window, 
not in the size of the window; the median of a window containing a single 
value is that single value; windows which contain no values are reported to 
have a median of 0.  For quantities like coverage which have a strong 
near-distance correlation within the input, this policy should be fine.


### Options

As this is an `awk` script, options may easily be set using a 
`parameter=value` syntax on the command line, for example:

```bash
windowWig windowsize=100 header=0 < input.dat | ...
```

Shown are commonly used parameters and their default values; check the `BEGIN`
section of the script for additional values that may be changed.

`header=1` : the number of header lines in the input (0 for none)

`skip_comment=1` : whether commend lines (beginning with '`#`') should be skipped (0 for no)

`no_data_value=0` : what should be printed for windows which have no data?

`windowsize=50` : the size of windows across which median values are calculated


### Caveats

Positions (column 2) must be sorted in increasing order within each 
reference (column 1).  They need not be consecutive. Note that the positions 
within each reference are assumed to be 
monotonically increasing in steps of 1 starting from 1 (by default) through the
last reported position within the reference, regardless of whether the data stream
actually contains values for every position.  So defined, every position from the
start to the end of every reference is guaranteed to be covered by a single 
reported window.

This script uses no `gawk` extensions.


[WIG]:  http://genome.ucsc.edu/goldenPath/help/wiggle.html



intervalBed
-----------

Starting with the same basic 3-column input as above, instead of continuous values
we have boolean 0/1 values in the third column (the header line is optional):

    REF         POS VAL
    reference1  1   0
    reference1  2   1
    reference1  3   1
    reference1  4   1
    reference1  5   1
    reference1  6   1
    reference1  7   1
    reference1  8   1
    reference1  9   0
    reference1  10  0
    reference1  11  1
    reference1  12  1
    reference1  13  1
    reference1  14  0
    ...

This script creates a [BED][] file defining intervals in which the value is true, note
that BED intervals are 0-based and [begin, end).  Output for the above is

    track name=intervalBed description="intervals of positions with 1s"
    reference1	1	8
    reference1	10	13

If the option `val_col=0` is given, then no third column of boolean values is
required, the presence of each position is given an implied 'true'.

The script also optionally allows for a grace distance through which the boolean
value need not be true for the interval to be maintained.  The grace distance will
only connect intervals, it will never begin or terminate a reference.  The output
from the above data with grace distance 10 is

    track name=intervalBed description="intervals of 1s, grace distance 10"
    reference1	1	13

Options can be changed with `option=value` on the command line.  The full set of
options available is

    FS, OFS      : input and output field separator characters (default "\t")
    header       : number of header lines on input (default 1)
    skip_comment : skip comment lines beginning with '#' (default 1, 'true')
    grace        : intervals separated by this much will be merged (default 0)
    min_width    : minimum width of an interval to produce (default 1)
    track        : print an initial track line on output? (default 1, 'true')
    trackname    : the value of 'name=' on the track line, default 'intervalBed'
    trackdesc    : the value of 'description=' on the track line, default set
                   from option values


Say you want to produce a track that noted the appearance of reads that mapped to
your reference at high quality.  This pipeline will pull out some high-quality 
pileup, produce a position-specific report of read mapping quality, keep only the 
reference, position and high-quality mapped-reads counts, turn the counts into 0/1 
using `boolify`, then uses this script to produce the bed track with a grace of
50bp:

```bash
samtools mpileup -sAB -d1000 -q20 -f ref.fa your.bam \
| smorgas --mapping-quality \
| cut -f1,2,5 \
| boolify col=3 header=1 \
| intervalBed header=1 grace=50 trackname=highQualityMappings \
      trackdesc="At least 1 high-quality-mapped read, grace 50bp"  > highqual.bed
```

Or a BED which marks intervals in which at least 10% of coverage is coming from
multiply-mapped (mapping quality 0) reads:

```bash
samtools mpileup -sAB -d1000 -q0 -f ref.fa your.bam \
| smorgas --mapping-quality \
| awk FS="\t" OFS="\t" '{ if (($4 / $3) >= 0.1) print; }' \
| cut -f1,2,4 \
| boolify col=3 header=1 \
| intervalBed header=1 trackname=multipleMappings \
      trackdesc="At least 10% multiply-mapped reads"  > dupmapped-10-percent.bed
```

[`smorgas`](https://github.com/douglasgscofield/smorgas) is under active development, and [`boolify`](https://github.com/douglasgscofield/tinyutils) is a super-simple script that turns the values in a column into 0 or 1.  Both are available here in their respective repositories.


### Options

As this is an `awk` script, options may easily be set using a 
`parameter=value` syntax on the command line, for example:

```bash
intervalBed grace=50 min_width=50 header=0 < input.dat | ...
```

Shown are commonly used parameters and their default values; check the `BEGIN`
section of the script for additional values that may be changed.

`header=1` : the number of header lines in the input (0 for none)

`skip_comment=1` : whether commend lines (beginning with '`#`') should be skipped (0 for no)

`grace=0` : set a grace distance (0 for none)

`min_width=1` : minimum interval width to be output

`track=1` : whether a BED-format `track` line should begin the output (0 for no)

`trackname="booleanIntervals"` : the name used in the `track` line

`trackdesc="intervals of 1s"` : the description used in the `track` line

`min_width=1` : whether 


### Caveats

Positions (column 2) must be sorted in increasing order within each 
reference (column 1).  They need not be consecutive. Note that the positions 
within each reference are assumed to be  monotonically increasing in steps of 1 
starting from 1 (by default) through the last reported position within the 
reference, regardless of whether the data stream actually contains values for 
every position.

This script uses no `gawk` extensions.


[BED]:  http://genome.ucsc.edu/FAQ/FAQformat.html#format1


samHeader2Bed.pl
----------------

Read a SAM header and produce a BED file or files from the reference sequence
descriptions.  The BED file contains reference sequences that satisfy filtering
criteria, such as minimum and/or maximum length, total length of reference
sequence included in each BED file, etc.  Options may be combined.

Say you want to call SNPs in contigs &ge; 1 kbp:

```bash
samtools view -H your.bam | samHeader2Bed.pl --min-length 1000 - > min1000.bed
samtools mpileup -l min1000.bed -u -f ref.fa your.bam | bcftools view ...
```

Or you want to gather pileups for successive ~10 Mbp chunks of contigs:

```bash
samtools view -H your.bam | samHeader2Bed.pl -o chunk --chunk-size 10000000 -
for BED in chunk.*.bed ; do
   samtools mpileup -l $BED -u -f ref.fa your.bam | gzip -c > $BED.mpileup.gz
done
```

**OPTIONS**

    -                     read SAM header from stdin, write to stdout
    -i FILE, --in FILE    read SAM header from FILE, else from stdin
    -o FILE, --out FILE   write BED output to FILE, else to stdout
    --num INT             include just the first INT sequences
    --chunk-size INT      produce separate BED files each containing sequences
                          representing approximately INT bp; reference sequences
                          are only described complete, so each BED is likely to 
                          describe more than INT bp.  Individual BED files are
                          named FILE.xx.bed where FILE is specified with --out, 
                          which is required, and xx is the integer sequence of 
                          file creation
    --min-length INT      do not include reference sequences shorter than INT
    --max-length INT      do not include reference sequences longer than INT
    --no-header           do not print header on output

    -?, --help            help message



fai2Bed.pl
----------

Identical in operation to `samHeader2Bed.pl` immediately above, except instead
of a SAM header the script creates a BED file from a Fasta index file in .fai
format.  Such a file is produced with `samtools faidx file.fa`, in this case
it would be named `file.fa.fai`.

All options are identical to `samHeader2Bed.pl`.



pileup2pro.pl
-------------

Convert pileup to profile format as used for input to mlRho.

**The developers of mlRho now provide a C-language tool `sam2pro` to
do the same job, available at <http://guanine.evolbio.mpg.de/mlRho/>.**

In profile format, the bases observed at each position in a mapping to a
reference sequence are enumerated.  This script simply converts the format, so
any filtering on base/mapping quality, etc. that you may wish to do should be
done when generating the pileup.

```bash
samtools mpileup -B -q1 -f ref.fa your.bam | pileup2pro.pl > mlRho-input.txt
```

mlRho (<http://guanine.evolbio.mpg.de/mlRho>) estimates population genetic
parameters from NGS data sequenced from a diploid genome.  See references below.

Profile format contains reference names, coordinates and raw numbers of bases:

    >contig_1
    1	0	2	0	0
    2	2	0	0	0
    3	0	2	0	0
    4	2	0	0	0
    5	0	0	0	2
    6	0	0	2	0
    7	2	0	0	0
    8	0	0	0	2
    9	0	0	0	2
    10	2	0	0	0
    ...

**OPTIONS**


    -                          read input from stdin, write to stdout
    --in FILE                  read input from FILE, else from stdin
    --out FILE                 write output to FILE, else to stdout
    --which-bams INT[,INT...]  produce profile output for the INT-th BAM(s) in 
                               order as provided on the samtools mpileup command
                               line, starting with 1; otherwise produce profile
                               output for all BAMs 
    --has-mapping-quality      must be specified if -s used for samtools mpileup
    --quiet                    don't print progress to stderr
    --help, -?                 help message

Pileup format created from multiple BAM files has 3 (or 4, with `-s`) columns per 
BAM file; this script will merge all columns while creating profile output up line
unless the `--which-bams` option is given.

Haubold B, P Pfaffelhuber, and M Lynch. 2010. mlRho - a program for estimating
the population mutation and recombination rates from sequenced diploid genomes.
*Molecular Ecology* 19s1:277-284.

Lynch M. 2008.  Estimation of nucleotide diversity, disequilibrium
coefficients, and mutation rates from high-coverage genome-sequencing projects.
*Molecular Biology and Evolution* 25:2421-2431.



mergePileupColumns
------------------

Merge columns of `samtools mpileup` listing pileup for several BAMs into a
single set of base call, base quality and (optionally) mapping quality columns.

**Note:** merging columns from separate `samtools mpileup` runs is not addressed by this script.

If multiple BAMs are given to a `mpileup` command, then the output includes 
separate coverage, base call and base quality scores for each BAM file.
With the `-s` option to `samtools mpileup`, the output for each BAM also 
includes a mapping quality score.  It is simple but not entirely straightforward 
to merge these columns into a single set of pileup columns.  Consider the 
following pileup line:

    ref1  1  A  4  ..,,  ffif  0  *  *  3  ,CC  gcd

The second BAM provides no coverage at this position, but if the `*` characters
are simply concatenated onto the base call and base quality strings, these two
columns now appear to indicate a gap and a valid base quality (with Phred+33 scaling)
or invalid base quality (with Phred+66 scaling), respectively, with counts that
no longer match the coverage:

    ref1  1  A  7  ..,,*,CC  ffif*gcd     <== incorrect merge
    ref1  1  A  7  ..,,,CC   ffifgcd      <== correct merge

Before samtools 1.0, an additional wrinkle wass introduced by an odd quirk in
`samtools mpileup -s` output.  If a BAM provides no coverage at a position, its
output does not include the mapping quality column.  Merging output from
`samtools mpileup -s`:

    ref1  1  A  4  ..,,  ffif  ]]]2  0  *  *  3  ,CC  gcd  +FF    <== missing column
    
    ref1  1  A  7  ..,,,CC   ffifgcd  ]]]2+FF                     <== correct merge

Starting with samtools 1.0, a fourth asterisk is output as expected.  The
script checks for '*' in a fourth column and skips four columns if it is there
(samtools 1.0+), and skips three columns if it is not (samtools < 1.0).

Yet another wrinkle occurs if there is coverage at this site, but it has all
been removed as a result of options to `samtools mpileup`, such as `--rf` to
remove some reads based on SAM flag values.  In these cases, coverage is 0, but
there are *four* columns of rather confused-looking output with empty base and
quality columns and quality characters in the mapping quality column.  This
apparently will be fixed at some point
(<http://seqanswers.com/forums/showthread.php?t=61759>).


### Usage

```bash
samtools mpileup -f ref.fa y1.bam y2.bam | mergePileupColumns > merged.pile
```

If the `-s` option was used for `samtools mpileup`, then set the `mpileup_s` variable
on the command line:


```bash
samtools mpileup -s -f ref.fa y1.bam y2.bam | mergePileupColumns mpileup_s=1 > merged.pile
```



extractFasta
------------

Extract named FASTA sequences quickly after creating a BioPerl FASTA database

SYNOPSIS

    extractFasta -n seqname2 -n seqname2 -i fasta.fa -o subset.fa
    extractFasta -fn file-of-seqnames -i fasta.fa -o subset.fa
    extractFasta -fn file-of-seqnames -i fasta.fa > subset.fa

Names must be given one per line in the names file.  Names of sequences
in the FASTA file are any characters after the initial '>' character in
the header of a FASTA sequence, followed by whitespace or an end of line,
and must be matched in their entirety.

OPTIONS

    -n NAME      name of FASTA sequence to extract (may be used multiple times)
    -fn FILE     file of FASTA sequence names, one per line
    -i FILE      input FASTA sequences
    -o FILE      output FASTA sequences (written to STDOUT if not specified)
    -I           do not remove intermediate index file created
    -?, --help   help message



extractFasta.pl
---------------

Extract a sequence or subsequence from a Blast database

USAGE: ./extractFasta.pl --db database --entry fasta-sequence-name [ --range L-R ]

First build the blast database using

    makeblastdb -parse_seqids -in sequence.fa -dbtype nucl|prot

with -dbtype dependent upon your input data.

Then you can extract the (sub)sequence of interest.  This is printed
to stdout in Fasta format.  If a subsequence is requested, the range
of the subsequence is added to the sequence identifier in the output.

OPTIONS:

    --db blast-db       The same name you would give to blastdbcmd,
                        which is exactly what this script does

    --entry sequence    Name of the sequence to extract from blast-db

    --range LOW-HIGH    1-based positions of a subsequence to extract
                        from the sequence, once it is found.




extractFastaSeqs.pl
-------------------

Extract named FASTA sequences.  Names must be given one per line in the names
file.  Names of sequences in the FASTA file are any characters after the
initial '>' character in the header of a FASTA sequence, followed by whitespace
or an end of line, and must be matched in their entirety.  Use `--header` to
match against the entire header of the FASTA sequence.  Use `--reverse` to
extract sequences that *do not match* any of the names.


All of these usages are equivalent:

```bash
extractFastaSeqs.pl --names subset-names.txt --in full.fa --out subset.fa
extractFastaSeqs.pl -n subset-names.txt -i full.fa - > subset.fa
cat full.fa | extractFastaSeqs.pl -n subset-names.txt - > subset.fa
extractFastaSeqs.pl subset-names.txt full.fa subset.fa
```

**OPTIONS**

    -                       read FASTA sequences from stdin and/or write
                            extracted sequences to stdout
    -n FILE, --names FILE   file containing names of FASTA sequences to extract
    -i FILE, --in FILE      input FASTA sequences
    -o FILE, --out FILE     output FASTA sequences
    --header                match entire contents of the FASTA header
    --reverse               output FASTA sequences that *do not match* any of
                            names given.
    --quit-on-seen          quit once all sequences in the names file are seen


    -?, --help              help message



trimFastq.pl
------------

Hard-trim a given number of bases from the 5' or 3' end (or both) from each read
in a file of FastQ-format reads, optionallly trimming to a maximum length.  Both
the sequence and quality strings are trimmed (naturally).

**OPTIONS**

    -                       read FastQ sequences from stdin and write FastQ
                            sequences to stdout
    --trim5 INT             trim INT bases from the 5' end of each read
    --trim3 INT             trim INT bases from the 3' end of each read
    --trimlen INT           trim from the 3' end so maximum length of each
                            read is INT


```bash
shuffleFastQ.pl r1.fq.gz r2.fq.gz - | trimFastq.pl --trim5 10 --trimlen 80 - | deshuffleFastQ.pl --minlen 50 - r1.trimmed.fq.gz r2.trimmed.fq.gz
```



shuffleFastq.pl, deshuffleFastq.pl
----------------------------------

Convert FastQ-format files from separate read 1 and read 2 files to interleaved
files, and back.  Input is compressed/decompressed automatically if the
filenames terminate with `.gz` or `.bz2`.  With the '`-`' option,
`shuffleFastq.pl` will write uncompressed FastQ to `/dev/stdout` and
`deshuffleFastq.pl` will read uncompressed FastQ from `/dev/stdin`.  With the
'`--md5`' option, MD5 checksums are also generated for the uncompressed FastQ,
see below for a further requirement for this to work.

```bash
shuffleFastq.pl  FA.1.fq FA.2.fq FA.i.fq.gz
deshuffleFastq.pl  FB.i.fq.gz FB.1.fq.gz FB.2.fq.gz
cat FC.i.fq | deshuffleFastq.pl - FC.1.fq.bz2 FC.2.fq.bz2
```

`shuffleFastq.pl` has a `--subset` option for specifying the maximum number of reads
to take from the beginning of the files.

`deshuffleFastq.pl` has a couple other options for filtering reads based
on minimum read length:

```bash
deshuffleFastq.pl --minlen 30 --single FD.se.fq.gz FD.i.fq.gz FD.1.fq.gz FD.2.fq.gz
```

This will only include read pairs in the output where both reads are at least
30bp long.  Any read that meets this criterion but has a mate that does not is
written to `FD.se.fq.gz`.  If the `--single` option is not specified, such reads are
dropped along with their mates.

Each also has a `--md5` option that will calculate MD5 checksums for the
uncompressed FastQ stream.  For each output file `outfile` (any `.gz` or `.bz2`
extension is removed), the checksum is saved in `md5sum` format to
`outfile.md5`; the filename in the checksum file is `outfile` without the
compression extension.  This requires the `compress_md5.sh` script from this
repository to be in your `PATH`.

These were originally built on the `shuffleSequences_fastq.pl` and
`deshuffleSequences_fastq.pl` scripts distributed with
[velvet](http://www.ebi.ac.uk/~zerbino/velvet) but have been completely
rewritten.



compress_md5.sh
---------------

Usage:

	compress_md5.sh  gz|bz2|none       filename
	compress_md5.sh  gz_self|bz2_self  filename
	
Saves stdin to filename (optionally compressed) while simultaneously computing
MD5 checksum on the uncompressed content, saved to `filename.md5`.

First argument is the compression format: `gz bz2 none`

If the compression format is `gz_self` or `bz2_self`, use filename as input
instead of stdin, and remove filename if compression is successful.

Second argument is the output filename without any compression suffix


Output is two files:
	
`filename`:
compressed if requested, with suffix as appropriate

`filename.md5`:
md5sum-format file with MD5 checksum for uncompressed content of the file; filename present in filename.md5 is filename without any compression suffix


checkFastqCollisions.pl
-----------------------

Usage:

    find . -type f -name '*.fastq.gz' | checkFastqCollisions.pl [ -v ] [ -t INT ]

Checks the list of filenames given on STDIN for collisions involving FastQ filenames.  If a collision is detected, a message is printed and the exit code of the script is non-zero.  The `find` tool is not required to have produced the filenames, but using it will avoid inadvertant duplicates on input.

For filenames ending with any combination of `.fastq` or `.fq` followed by `.gz` or `.bz2` or `.xz`, three types of collisions are detected:

* Identical filenames after all directory information is removed
* Identical directory/filenames when keeping the final directory before the filename
* Identical complete filenames specified on input, which should be avoided as it will create both of the above identities

With the `-v` or `--verbose` option, the complete filename of each duplicate is also printed.

With the `-t INT` or `--trimsuffix INT` option, then collisions are also detected after `INT` dot-suffixes are removed from the filename.  Using this option, a collision between `file.fastq.gz` and `file.fastq.bz2` is detected with `-t 1`, and between `file.fastq.gz` and `file.fq.bz2` with `-t 2`.


fastaOneline
------------

Usage:

    fastaOneline  [ -1 ] file.fa  >  output.fa

Prints all sequences from file.fa as single lines.  The output is Fasta
format (unless -1 is used), but there is no wrapping of sequence lines.
This can be useful for grep and count operations.  Uses BioPerl.

**Options**

    -1     Truly a single line: sequence-name <TAB> sequence <NEWLINE>


fastaGC.pl
----------

Usage:

    fastaGC.pl [options] [file ...]

By default `fastaGC.pl` will print length, GC and base content information
across all FASTA sequences in the input file(s). More detailed
information may be requested with options. For all information, GC
content is printed in the form of a tab-separated table suitable for
import into whatever other program can read tables.

Some options can print GC content on a block-by-block basis. These
options break input sequences into blocks, where each block is one line
of FASTA sequence. The scale of a block for these options is dependent
upon the input format. If sequences are broken into 60-character lines,
then block size will be 60 bp.

The options are not mutually exclusive. All four types of GC content may
be displayed in a single run.

One advantage of this script is that it doesn't use BioPerl.  Some
disadvantages are that it is slow, could use a streamlining, and it would be
nice if it produced full-composition output similar to Jim Kent's `faCount`
program.

**Options**

    --nototal    do not print information for total input
    --concat     print block-by-block GC content for input sequences
                 as if they were one concatenated sequence
    --block      print mean GC content across all input sequences
                 on a block-by-block basis
    --seq        print GC content for each individual sequence
                 in the input file(s)
    --verbose    print progress (every 1000 sequences), and print
                 headers prior to the output of each table of GC content
    --help, -?   brief help message




gmhmmp2Table.pl
---------------

    gmhmmp2Table.pl file.in > file.out


Read the output of [gmhmmp](http://www.genepro.com/Manuals/EuGM/EuGM_usage.aspx), an ORF-prediction program, and produce a table summarizing each of the predicted ORFs.


gmhmmp2Fasta.pl
---------------

    gmhmmp2Fasta.pl file.in > file.out


Read the output of [gmhmmp](http://www.genepro.com/Manuals/EuGM/EuGM_usage.aspx), an ORF-prediction program, and produce a Fasta file containing the sequence of each of the predicted ORFs.


convertSequence.pl
------------------

Sequence format conversion - format to format - all available in BioPerl's [Bio::SeqIO](http://bioperl.org/howtos/SeqIO_HOWTO.html)

Usage: `convertSequence.pl [ --if input-format ] [ [--input] input.file ] [ --of output-format ] [ --output output.file ]`

    --input     infile    default from STDIN, or first non-option command line argument
    --output    outfile   default to STDOUT

    --if|format fasta     input format, see Bio::SeqIO documentation for available formats
                genbank   http://bioperl.org/howtos/SeqIO_HOWTO.html
                swiss
                embl      default genbank
                etc.

    --of|ormat  fasta     output file format, as above for formats
                etc.      default fasta

    -v|--verbose          print end summary of conversion effort
    -h|--help             this help


convertAlignment.pl
------------------

Alignment format conversion - format to format - all available in BioPerl's [Bio::AlignIO](http://bioperl.org/howtos/AlignIO_and_SimpleAlign_HOWTO.html).

Can allow degap the output with `--degap`.

Usage: `convertAlignment.pl [ --if input-format ] [ [--input] input.file ] [ --of output-format ] [ --output output.file ]`

    --input     infile    default from STDIN, or first non-option command-line argument
    --output    outfile   default to STDOUT

    --if|format fasta     input format, see Bio::AlignIO documentation for available formats
                clustalw  http://bioperl.org/howtos/AlignIO_and_SimpleAlign_HOWTO.html
                nexus
                phylip    default clustalw
                etc.

    --of|ormat  fasta     output file format, as above for formats
                etc.      default fasta

    --degap               convert aligned sequences to standard sequences by removing gaps
    --uc | --uppercase    convert sequence to uppercase


    -v | --verbose        print end summary of conversion effort
    -h | --help           this help


cutadaptReportScript.sh
-----------------------

Usage:

    cutadaptReportScript.sh reads-a.cutReport [ reads-b.cutReport ... ] > cutadapt_report_table.txt


Collect the results of `*.cutReport` files produced by the [cutadapt](https://code.google.com/p/cutadapt/) adapter-trimming tool.  This is a useful way to quickly assess which are the set of adapters used in different read sets.  Running `cutadapt` requires a set of putative adapter sequences, and it is the results for that specific set that are reported with this script.  Written in collaboration with Amaryllis Vidalis in the EMG department at Umeå University.  **Requires the same set of adapters are specified in the same order for each collected run.**


coords-view
-----------

Usage:

    nucmer ... # creates out.delta
    show-coords -lcdTH out.delta > out.sc
    coords-view out.sc

Script for creating PNG visualising the results of a nucmer/mummer alignment.
From Alan Twaddle, Sourceforge repository at
<http://sourceforge.net/projects/coordsview/>.  Modified a little bit by me.


stacksSNPAssays.pl
------------------

Mine Stacks output to design sequence templates to use when designing SNP
assays (iPLEX, iSelect).  Configurable.  Under construction.



pos2bed
-------

Convert chr,pos file to BED intervals.

    --posfile FILE   Two-column file with the first two columns specifying chromosome,position locations.
                     Position is 1-based, so the posfile is like the first two columns of a VCF or GFF/GTF
                     file.  Only the first two columns are used, if other columns are present they are ignored.
                     Lines beginning with '#' are ignored.
    --faifile FILE   FAI-format file describing sizes of chromosomes expected in the posfile.  Lines beginning
                     with '#' are ignored.  If --no-addends is specified, this file is not required.  
    --no-addends     Do not add the ends of chromosomes.  By default, BED intervals are created for the
                     beginning and end of each chromosome, with the beginning at position 1 and the end at the
                     length given in the faifile.  When this option is specified, these intervals are not 
                     created, and the --faifile option is not required.
    --help           This help

