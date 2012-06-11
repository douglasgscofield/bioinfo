pileVar.pl
----------

Call variants by processing output of samtools mpileup command.

The initial development of this script is for the calling of variants during
reference-guided assembly of a chloroplast sequence.  Thus, my first interest
is uncovering *fixed differences* between the reference and the set of reads,
and identifying regions of ambiguity due to indels, rather than finding SNPs
isolated against a background of sequence variation.  Eventually I may extend
this or write another script to handle SNP calling in pooled samples.

Built upon the structure of the very useful mpileup-format parser script
available from Galaxy at
<https://bitbucket.org/galaxy/galaxy-central/src/tip/tools/samtools/pileup_parser.pl>.

### Input

Input is **raw** mpileup output, via something like

    samtools mpileup -BQ0 -d 1000000 -f ref.fa -D your.bam -A > mpileup.raw.txt

This requests no read filtering, and includes anomalous base pairs (-A).  You
might want to drop that last one or include others, based on your pipeline.


### Output

Current output is a table summarizing base counts and qualities, annotated with
notes about indels and likely fixed positions.  Various thresholds and other
aspects of output are controlled with a number of options.

    samtools mpileup -BQ0 -d 1000000 -f ref.fa -D your.bam -A | ./pileVar.pl

will produce

    sequence	coord	ref	cvg	A	C	G	T	qualA	qualC	qualG	qualT	qual_cvg	ref_diff	cons	cons_qual	note
    T_CP|NC_009143.1	a	1	10	9	0	0	1	35.0	0	0	36.0	10	1	A	35.0
    T_CP|NC_009143.1	t	2	21	0	0	0	21	0	0	0	36.2	21	0	T	36.2
    T_CP|NC_009143.1	g	3	41	0	0	41	0	0	0	36.6	0	41	0	G	36.6
    T_CP|NC_009143.1	g	4	49	0	0	49	0	0	0	36.7	0	49	0	G	36.7
    T_CP|NC_009143.1	g	5	58	0	0	58	0	0	0	35.9	0	58	0	G	35.9
    ...
    T_CP|NC_009143.1	t	56	729	0	0	0	729	0	0	0	36.7	729	0	T	36.7
    T_CP|NC_009143.1	g	57	748	0	0	748	0	0	0	36.7	0	748	0	G	36.7
    T_CP|NC_009143.1	a	58	739	739	0	0	0	36.3	0	0	0	739	0	A	36.3
    T_CP|NC_009143.1	t	59	751	0	0	653	98	0	0	36.4	35.8	751	653	G	36.4	fixed T -> G (0.870)
    T_CP|NC_009143.1	c	60	757	0	757	0	0	0	36.9	0	0	757	0	C	36.9
    T_CP|NC_009143.1	c	61	767	0	767	0	0	0	36.9	0	0	767	0	C	36.9
    T_CP|NC_009143.1	a	62	779	777	0	0	2	35.3	0	0	39.0	779	2	A	35.3
    ...
    T_CP|NC_009143.1	a	180	1257	1250	2	3	1	37.4	34.5	38.3	39.0	1256	6	A	37.4
    T_CP|NC_009143.1	a	181	1304	1275	6	3	20	37.5	37.4	36.0	35.5	1304	29	A	37.5
    T_CP|NC_009143.1	a	182	1472	1459	2	3	8	36.7	36.0	39.0	35.6	1472	13	A	36.7
    T_CP|NC_009143.1	g	183	1541	305	16	1182	36	36.4	36.4	37.0	36.5	1539	357	G	37.0
    T_CP|NC_009143.1	t	184	1578	366	16	25	1170	37.1	36.7	35.7	36.0	1577	407	T	36.0	indel fraction +0.061 -0.001
    indel, +, 1 bases, +1A, 4 times
    indel, +, 1 bases, +1C, 1 times
    indel, +, 1 bases, +1G, 3 times
    indel, +, 1 bases, +1g, 1 times
    indel, +, 2 bases, +2AA, 13 times
    indel, +, 2 bases, +2aa, 7 times
    indel, +, 3 bases, +3AAA, 42 times
    indel, +, 3 bases, +3aaa, 26 times
    indel, -, 1 bases, -1A, 1 times
    T_CP|NC_009143.1	a	185	1649	1628	3	6	11	37.5	36.4	36.6	36.6	1648	20	A	37.5
    T_CP|NC_009143.1	a	186	1677	1668	4	2	3	37.5	38.3	37.1	33.1	1677	9	A	37.5
    ...


### Options

Output of `./pileVar.pl --help`:

~~~~

NAME

  pileVar.pl - parse SAMtools mpileup format for variant information

SYNOPSIS

  pileVar.pl  [options]

OPTIONS

    -                        read input from STDIN, write to STDOUT
    -in <filename>           read input from <filename>, else from STDIN
    -out <filename>          write output to <filename>, else to STDOUT
    -noheader                do not print header on output
    -base_qual_offset <int>  offset of base quality ASCII value from 0 quality [default 64]

    -ploidy <int>            ploidy of samples from which pileup is derived [default 1]

    -no_coord_check          do not check coordinates for increase-by-1 consistency.  If a consensus
                             FASTA file is being produced, this suppresses the insertion of gaps
                             in the output
    -variants_only           print out only those coordinates containing variants from the reference
                             [default 0]
    -consensus               call consensus sequence from reads, adds two columns to output,
                             one for consensus call, one for mean quality of call [default 1]
    -consensus_fasta <filename>  print the consensus sequence in FASTA format to the given file
    -fasta_gap_char <char>   the character used for coordinate-skipping gaps in the consensus
                             FASTA file, gaps are also identified by position in the FASTA header line 
                             [default n]
    -print_bases_quals       print out bases and base qualities from input mpileup
                             [default 0]

  Quality:

    -mincov <int>            minimum raw coverage to call a change to reference [default 0]
    -minqual <int>           minimum base quality cutoff [default 0]
    -minqualcov <int>        minimum coverage of quality bases to call a change to reference [default 0]
    -qual                    print qualities of variants [default 1]
    -qualcalc mean|rms|sum   method for calculating qualities of variants [default rms]
    -qualround <int>         digit to which to round variant qualities [default 1]

  Indels:

    -track_indels            track the presence of indels [default 0]
    -filter_indel_frac <fraction>  do not report indels at a coordinate if the fraction of 
                             reads containing them is below <fraction> [default 0.05]

  SNP variants:

    -fixed_frac <fraction>   minimum fraction of non-ref bases required to call a base as 
                             fixed with respect to the reference [default 0.8]

    -help, -?                help message

~~~~

