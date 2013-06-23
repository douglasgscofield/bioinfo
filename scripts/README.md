Bioinformatics scripts
======================

These have been useful to me, I hope they can be useful to you!

mummer2Vcf.pl
-------------

Convert a list of SNPs and indels produced by [Mummer's](http://mummer.sourceforge.net/)
`show-snps -T` command to a rough substitute for a 
[VCF file](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41).
VCF requires the first reference base to be output for indels and Mummer
doesn't output this, so this script requires a file containing the Fasta-format
reference sequence, in SNP order.  It will open the file named as the reference
in the first line of the `show-snps -T`-formatted file, but this can also be
supplied with the `--fasta` option, which will override the file-encoded
reference sequence if it is provided.

It also provides the `--type SNP|INDEL` option to subselect the set of variants
produced if you like, as well as the `--snpEffect` option to produce variants suitably formatted
for the still-in-infancy [`snpEffect()`](https://github.com/douglasgscofield/snpEffect) function.

The 'VCF' format produced by this script is not actually VCF yet, but I'm working
on it :-)

This script requires BioPerl 1.6.1.


phredDetector.pl
----------------

Heuristically determine the [Phred-scaled quality score](http://en.wikipedia.org/wiki/FASTQ_format) used in the FastQ file presented on input.  If everything looks reasonable, it simply prints `33`, `64` or `59` (this last for older Solexa sequences) to `stdout`, so it can be used in pipeline scripts to autodetect Phred encodings:

```bash
FastQ_file="fastq_file.fq.gz"
Quality=$(phredDetector.pl $FastQ_file)
if [ "$Quality" = "33" ] ; then
   echo "perhaps Illumina 1.8+"
elif [ "$Quality" = "64" ] ; then
   echo "perhaps Illumina 1.3+ to pre-1.8"
elif [ "$Quality" = "59" ] ; then
   echo "perhaps Solexa"
else
   echo "Couldn't autodetect quality for $FastQ_file, return value was '$Quality'"
fi
```


BioPerl is not required.  Usage is:

    phredDetector.pl  [ options ]  fastq_file1.fq[.gz]

Input must be in FastQ format, if a filename is given it may be gzipped (*.gz)

Output to stdout is either

* `33`, appears to be Illumina 1.8+ or Sanger quality encoding
* `64`, appears to be Illumina 1.3+ to pre-1.8 quality coding
* `59`, appears to be Solexa (pre-Illumina) base-64 quality

The only data interpreted are read quality scores on line 4 of each read; all sequence, pairing information etc. is ignored.

**Caveat:** If all bases on input are of unusually high quality, then a Phred base of 59 or
64 may be reported when a Phred scale based on 33 was the one actually used.  A
few heuristics are used to detect possible problems, but these are not
comprehensive.

* If (maximum quality &minus; minimum quality) &ge; 47, a warning message is printed
  to stderr and detected quality encoding (possibly erroneous) to stdout; this value can be adjusted with the `--wide` option
* If (maximum quality &minus; minimum quality) &le; 20, a warning message is printed
  to stderr and detected quality encoding (possibly erroneous) to stdout; this value can be adjusted with the `--narrow` option
* If otherwise unusual quality scores or unknown input were detected, an error
  message is printed to stderr and '??' to stdout

This script does not diagnose faulty FastQ files, nor does it fully diagnose
the various versions of Solexa, Sanger, Illumina pipelines described in
http://en.wikipedia.org/wiki/FASTQ_format.  It simply applies a decreasing minimum cutoff
from 64, to 59, to 33.  It is thus compatible with Sanger encoding but
it does not take into account finer distinctions such as Illumina 1.3+ pipelines
not typically producing quality scores 0 and 1.

**Options**

    -             Read uncompressed FastQ from stdin
    --reads INT   Number of reads to process to determine Phred basis [10000]
                  If 0, process *all* reads in the input file
    --wide INT    Use INT for the 'too wide' first heuristic above [47]
    --narrow INT  Use INT for the 'too narrow' second heuristic above [20]

    --help | -?   Generate this help output


fermiExtractContigs.pl
------------------------

Extract Fasta-format contigs from a
[fermi](https://github.com/lh3/fermi)-format `*.fq.gz` FastQ-like scaftig
files.  Writes Fasta to `stdout`, with each sequence given its fermi sequence
name, and a description that includes sequence length, number of non-redundant
reads that built the scaftig, and median coverage of non-redundant reads along
the scafftig.  So this:

````
@26417937:25351227_0	191	83
TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTATATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTAAACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCCTGCATGTAGAA
+
"#$%%&'((()**+,-../0123456789:;<=>?@ABCDEFGHIJKLLMNOPQRSTUVWXYZ[\]^_`abcdeffghijklmnoopqrstttttttttttsrqpponmmmlkkjihggfedcba`_^]\[ZYXWVUTSRQPONMLKJIIHGFEDCBA@?>=<;:9876543210//.-,+*)('&&%$#"
````

Becomes this:

````
>26417937:25351227_0 length:191,n_reads:83,median_coverage:44
TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTAT
ATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTA
AACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCC
TGCATGTAGAA
````

This script required BioPerl 1.6.1.



fermiExtractContigs_simple.sh
------------------------

Extract Fasta-format contigs from a
[fermi](https://github.com/lh3/fermi)-format `*.fq.gz` FastQ-like scaftig
files.  Writes Fasta to `stdout`, with each sequence given its fermi sequence name and a description that includes sequence length and number of non-redundant reads that built the scaftig.  So this:

````
@26417937:25351227_0	191	83
TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTATATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTAAACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCCTGCATGTAGAA
+
"#$%%&'((()**+,-../0123456789:;<=>?@ABCDEFGHIJKLLMNOPQRSTUVWXYZ[\]^_`abcdeffghijklmnoopqrstttttttttttsrqpponmmmlkkjihggfedcba`_^]\[ZYXWVUTSRQPONMLKJIIHGFEDCBA@?>=<;:9876543210//.-,+*)('&&%$#"
````

Becomes this:

````
>26417937:25351227_0 length:191,n_reads:83
TTTCTATTCTAAACCACCGTATATATGTAATTTCTATTCTAAACTAACCTGTGTCCGTATATATGTAATTTCTATTCTAAACTACCTGTGTGAAGAAGCCCTACGTTTCTTTCTATTCTAAACTACCGTATTTCCTTACGTTTTTTTCTATTCTTTTCCACTCAAAATGGCCGACACTCCTGCATGTAGAA
````

**Differences from `fermiExtractContigs.pl`**: In the interests of speed and simplicity, this script does no wrapping of the fermi-generated sequences, likewise it does not include median coverage along the entire length of the sequence.  If these limitations are not suitable for your task, consider using the `fermiExtractContigs.pl` script above.


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

````bash
windowWig windowsize=100 header=0 < input.dat | ...
````

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

````bash
intervalBed grace=50 min_width=50 header=0 < input.dat | ...
````

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

````bash
samtools view -H your.bam | samHeader2Bed.pl --min-length 1000 - > min1000.bed
samtools mpileup -l min1000.bed -u -f ref.fa your.bam | bcftools view ...
````

Or you want to gather pileups for successive ~10 Mbp chunks of contigs:

````bash
samtools view -H your.bam | samHeader2Bed.pl -o chunk --chunk-size 10000000 -
for BED in chunk.*.bed ; do
   samtools mpileup -l $BED -u -f ref.fa your.bam | gzip -c > $BED.mpileup.gz
done
`````

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



pileup2pro.pl
-------------

Convert pileup to profile format as used for input to mlRho.  This is a format
in which bases present at each position in a mapping to a reference sequence
are enumerated.  This script simply converts the format, so any filtering on
base/mapping quality, etc. that you may wish to do should be done when
generating the pileup.

````bash
samtools mpileup -B -q1 -f ref.fa your.bam | pileup2pro.pl > mlRho-input.txt
````

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

An additional wrinkle is introduced by an odd quirk in `samtools mpileup -s` output.
If a BAM provides no coverage at a position, its output does not include the
mapping quality column.  Merging output from `samtools mpileup -s`:

    ref1  1  A  4  ..,,  ffif  ]]]2  0  *  *  3  ,CC  gcd  +FF    <== missing column
    
    ref1  1  A  7  ..,,,CC   ffifgcd  ]]]2+FF                     <== correct merge


### Usage

````bash
samtools mpileup -f ref.fa y1.bam y2.bam | mergePileupColumns > merged.pile
````

If the `-s` option was used for `samtools mpileup`, then set the `mpileup_s` variable
on the command line:


````bash
samtools mpileup -s -f ref.fa y1.bam y2.bam | mergePileupColumns mpileup_s=1 > merged.pile
````



extractFastaSeqs.pl
-------------------

Extract named FASTA sequences.  Names must be given one per line in the names
file.  Names of sequences in the FASTA file are any characters after the
initial '>' character in the header of a FASTA sequence, followed by whitespace
or an end of line, and must be matched in their entirety.  Use `--header` to
match against the entire header of the FASTA sequence.  Use `--reverse` to
extract sequences that *do not match* any of the names.


All of these usages are equivalent:

````bash
extractFastaSeqs.pl --names subset-names.txt --in full.fa --out subset.fa
extractFastaSeqs.pl -n subset-names.txt -i full.fa - > subset.fa
cat full.fa | extractFastaSeqs.pl -n subset-names.txt - > subset.fa
extractFastaSeqs.pl subset-names.txt full.fa subset.fa
`````

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



shuffleFastq.pl, deshuffleFastq.pl
----------------------------------

Convert FastQ-format files from separate read 1 and read 2 files to interleaved
files, and back.  Input is automatically un-gzipped if required with no
particular filename requirements, and output is gzipped if the filename ends
with `.gz`.  All the gzipping/ungzipping is done in pipes.  With the '`-`' option,
`shuffleFastq.pl` will write to `stdout` and `deshuffleFastq.pl` will read from
`stdin`.

````bash
shuffleFastq.pl  FA.1.fq FA.2.fq FA.i.fq.gz
deshuffleFastq.pl  FB.i.fq.gz FB.1.fq.gz FB.2.fq.gz
cat FC.i.fq | deshuffleFastq.pl - FC.1.fq.gz FC.2.fq.gz
````

etc.  `deshuffleFastq.pl` has a couple other options for filtering reads based
on minimum read length:

````bash
deshuffleFastq.pl -minlen 30 -single FD.se.fq.gz FD.i.fq.gz FD.1.fq.gz FD.2.fq.gz
````
  
This will only include read pairs in the output where both reads are at least
30bp long.  Any read that meets this criterion but has a mate that does not is
written to `FD.se.fq.gz`.  If the `-single` option is not specified, such reads are
dropped along with their mates.

These were originally built on the the `shuffleSequences_fastq.pl` and
`deshuffleSequences_fastq.pl` scripts distributed with
[velvet](http://www.ebi.ac.uk/~zerbino/velvet).


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

One advantage of this script is that it doesn't use BioPerl.  A disadvantage
of this script is that it could use a streamlining, and would be more
useful if it produced full-composition output similar to Jim Kent's `faCount`
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

Usage:

    convertSequence.pl -f genbank file.in -of fasta - > file.out


Read the sequence data in GenBank format from `file.in`, convert it to Fasta format and write it via `stdout` to `file.out`.


convertAlignment.pl
------------------

Usage:

    convertAlignment.pl -f clustalw file.in -of phylip - > file.out


Read the multiple sequence alignment in ClustalW format from `file.in`, convert it to [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) format and write it via `stdout` to `file.out`.


gmhmmp2Fasta.pl
---------------

Usage:

    gmhmmp2Fasta.pl file.in > file.out


Read the output of [gmhmmp](http://www.genepro.com/Manuals/EuGM/EuGM_usage.aspx), an ORF-prediction program, and produce a Fasta file containing the sequence of each of the predicted ORFs.


cutadaptReportScript.sh
-----------------------

Usage:

    cutadaptReportScript.sh reads-a.cutReport [ reads-b.cutReport ... ] > cutadapt_report_table.txt


Collect the results of `*.cutReport` files produced by the [cutadapt](https://code.google.com/p/cutadapt/) adapter-trimming tool.  This is a useful way to quickly assess which are the set of adapters used in different read sets.  Running `cutadapt` requires a set of putative adapter sequences, and it is the results for that specific set that are reported with this script.  Written in collaboration with Amaryllis Vidalis in the EMG department at Umeå University.  **Requires the same set of adapters are specified in the same order for each collected run.**
