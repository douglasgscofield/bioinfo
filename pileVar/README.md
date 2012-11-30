pileVar.pl
==========

Summarize variants in `samtools mpileup` output.

The initial development of this script was for the identification of variants 
during reference-guided assembly of a chloroplast sequence. So my first 
interest was to uncover *fixed differences* between the reference and the set of reads,
and identifying regions of ambiguity due to indels, in contrast to finding SNPs and
other variants isolated against a background of sequence variants.

It has since evolved into a bit more of a general-purpose script for summarizing
the raw "model-free" (or with options, "naive-model") variation present within pileup. 
Other tools do a good job of calling SNPs and indels (sort of) and this script 
could be thought of as providing a null-model summary of all variants present, to
allow for judging whether the variants it shows should or should not be called by
other SNP- and indel-calling tools.

Currently it prints its own simple reports, but it would be nice to have it print
VCF to allow for comparisons.

The pileup-parsing code and indeed the entire script was initially built upon the 
structure of the very useful mpileup-format parser script available from Galaxy at
<https://bitbucket.org/galaxy/galaxy-central/src/tip/tools/samtools/pileup_parser.pl>.
It's moving away from that as I replace more and more of it, but that script gave me
a great head start while I was learning about pileup.


Input
-----

Input is `samtools mpileup` output, via something like

````bash
samtools mpileup -AB -q0 -Q0 -d1000000 -f ref.fa -D your.bam | pileVar.pl ...
````

This is pretty close to raw pileup, without (all?) the quality adjustments
that `samtools` will do with BAQ computation etc. It's a good place to start
if you want to use this script as discussed above, to do a model-free summary
of variation within your mapped reads.


Output
------

#### --indel-mode

The output flavor I'm using currently is `--indel-mode`, which summarizes 
all indels present in the pileup and accompanies each summary with a good/bad 
tag based on very simple criteria:

* all indel lengths and indel operations are identical (the `ilen` and `ioper` 
  columns), after converting operations for reverse-orientation reads to 
  uppercase
* the indel frequency is greater than or equal to the `--indel-frac` frequency
  (below uses `--indel-frac 0.1`)

Example output:

    seq     pos     ref     cvg     iqual   ifreq   ilen    ioper   n_var
    seq1    1749    C       30      bad     0.0333  -4      -4GAAA  1
    seq1    1750    G       26      bad     0.0769  -1      -1A     2
    seq1    1758    A       25      bad     0.0400  3       +3TGG   1
    seq1    1771    G       32      bad     0.0312  -1      -1T     1
    seq1    1773    T       32      bad     0.0312  1       +1C     1
    seq1    1776    C       35      bad     0.0857  1,-1    +1G,-1T 3
    seq1    2076    A       54      bad     0.0741  -1      -1G     4
    seq1    2078    A       55      bad     0.0727  -1      -1G     4
    seq1    2281    A       29      bad     0.0345  1       +1T     1
    seq1    2294    C       32      bad     0.0312  4       +4TGGG  1
    seq1    3408    T       23      good    0.8261  -1      -1N     19
    seq1    4386    A       16      good    0.3750  -29     -29GACAGAAAGATGCTAAGGACGGATTTAGC        6
    seq1    5262    G       49      bad     0.0204  -1      -1A     1
    seq1    5369    A       49      bad     0.0204  1       +1G     1
    seq1    5393    G       55      bad     0.0182  1       +1A     1
    seq1    5821    G       237     bad     0.0042  -2      -2GA    1
    seq1    5849    G       167     good    0.1138  1       +1A     19
    ...

You could quickly generate a distribution of indel sizes with

````bash
samtools mpileup -AB -q1 -f ref.fa -D your.bam \
    | pileVar.pl --indel-mode --indel-frac 0.1 \
    | grep '\bgood\b' | cut -f7 \
    | hist
````

`hist` is a little awk script that generates a histogram, if I haven't posted 
it and you'd like it, [drop me a line](mailto:douglasgscofield@gmail.com).


#### Default output

Default output is a table summarizing base counts and qualities, annotated with
notes about indels and likely fixed positions. Various thresholds and other
aspects of output are controlled with a number of options.

````bash
samtools mpileup -BQ0 -d 1000000 -f ref.fa -D your.bam -A | ./pileVar.pl
````

will produce

    sequence	coord	ref	cvg	A	C	G	T	qualA	qualC	qualG	qualT	qual_cvg	ref_diff	cons	cons_qual	note
    ref.fa	a	1	10	9	0	0	1	35.0	0	0	36.0	10	1	A	35.0
    ref.fa	t	2	21	0	0	0	21	0	0	0	36.2	21	0	T	36.2
    ref.fa	g	3	41	0	0	41	0	0	0	36.6	0	41	0	G	36.6
    ref.fa	g	4	49	0	0	49	0	0	0	36.7	0	49	0	G	36.7
    ref.fa	g	5	58	0	0	58	0	0	0	35.9	0	58	0	G	35.9
    ...
    ref.fa	t	56	729	0	0	0	729	0	0	0	36.7	729	0	T	36.7
    ref.fa	g	57	748	0	0	748	0	0	0	36.7	0	748	0	G	36.7
    ref.fa	a	58	739	739	0	0	0	36.3	0	0	0	739	0	A	36.3
    ref.fa	t	59	751	0	0	653	98	0	0	36.4	35.8	751	653	G	36.4	fixed T -> G (0.870)
    ref.fa	c	60	757	0	757	0	0	0	36.9	0	0	757	0	C	36.9
    ref.fa	c	61	767	0	767	0	0	0	36.9	0	0	767	0	C	36.9
    ref.fa	a	62	779	777	0	0	2	35.3	0	0	39.0	779	2	A	35.3
    ...
    ref.fa	a	180	1257	1250	2	3	1	37.4	34.5	38.3	39.0	1256	6	A	37.4
    ref.fa	a	181	1304	1275	6	3	20	37.5	37.4	36.0	35.5	1304	29	A	37.5
    ref.fa	a	182	1472	1459	2	3	8	36.7	36.0	39.0	35.6	1472	13	A	36.7
    ref.fa	g	183	1541	305	16	1182	36	36.4	36.4	37.0	36.5	1539	357	G	37.0
    ref.fa	t	184	1578	366	16	25	1170	37.1	36.7	35.7	36.0	1577	407	T	36.0	indel fraction +0.061 -0.001
    indel, +, 1 bases, +1A, 4 times
    indel, +, 1 bases, +1C, 1 times
    indel, +, 1 bases, +1G, 3 times
    indel, +, 1 bases, +1g, 1 times
    indel, +, 2 bases, +2AA, 13 times
    indel, +, 2 bases, +2aa, 7 times
    indel, +, 3 bases, +3AAA, 42 times
    indel, +, 3 bases, +3aaa, 26 times
    indel, -, 1 bases, -1A, 1 times
    ref.fa	a	185	1649	1628	3	6	11	37.5	36.4	36.6	36.6	1648	20	A	37.5
    ref.fa	a	186	1677	1668	4	2	3	37.5	38.3	37.1	33.1	1677	9	A	37.5
    ...


### Options

Output of `pileVar.pl --help`:

````
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
    --base-qual-offset INT   offset of base quality ASCII value from 0 quality [default 64]

    --ignore-N               do not produce output for positions with reference of N

    --ploidy INT             ploidy of sample, currently ignored [default 2]

    --no-pos-check           do not check positions for increase-by-1 consistency.  If a consensus
                             FASTA file is being produced, this suppresses the insertion of gaps
                             in the output (--no-coord-check is a synonym)
    --hetz                   determine heterozygous positions, applying all other options
    --hetz-min-freq FLOAT    minimum frequency for heterozygous call at a site [0.1]
    --hetz-check             check that #alleles do not exceed ploidy, if so ignore it [1]
                             this all would be better replaced with a model
    --variants-only          print out only those positions containing variants from the reference
                             [default 0]
    --consensus              call consensus sequence from reads, adds two columns to output,
                             one for consensus call, one for mean quality of call [default 0]
    --consensus-fasta FILE   print the consensus sequence in FASTA format to the given file
    --fasta-gap-char CHAR    the character used for position-skipping gaps in the consensus
                             FASTA file, gaps are also identified by position in the FASTA header line 
                             [default n]
    --print-bases-quals      print out bases and base qualities from input mpileup
                             [default 0]

  Quality:

    --mincov INT             minimum raw coverage to call a change to reference [default 0]
    --minqual INT            minimum base quality cutoff [default 0]
    --minqualcov INT         minimum coverage of quality bases to call a change to reference [default 0]
    --qual                   print qualities of variants [default 0]
    --qualcalc mean|rms|sum  method for calculating qualities of variants [default rms]
    --qualround INT          digit to which to round variant qualities [default 1]

  Indels:

    --indels                 track the presence of indels [default 0]
    --indel-frac FLOAT       do not report indels at a position if the fraction of 
                             reads containing them is below FRAC [default 0.05]
    --indel-mode             track ONLY the presence of indels [default 0]

  SNP variants:

    --fixed-frac FLOAT       minimum fraction of non-ref bases required to call a base as 
                             fixed with respect to the reference [default 0.8]

    --help, -?               help message

````
