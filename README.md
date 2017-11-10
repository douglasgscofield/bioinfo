Bioinformatics tools
====================

More details are provided within each subfolder.

### [pileVar](https://github.com/douglasgscofield/bioinfo/tree/master/pileVar)

**pileVar.pl** identifies variants based on mpileup output from
[samtools](http://samtools.sourcefourge.net).  This tool is currently geared
toward identifying fixed differences and regions of ambiguity due to indels.


### [plotGFF](https://github.com/douglasgscofield/bioinfo/tree/master/plotGFF)

Plot simple schematic of annotation features described in a GFF file.


### [plotSpatialPies](https://github.com/douglasgscofield/bioinfo/tree/master/plotSpatialPies)

Plot spatially-distributed pie charts.


### [R-scripts](https://github.com/douglasgscofield/bioinfo/tree/master/R-scripts)

A handful of bioinformatics scripts in R.


### [stacksTemplates](https://github.com/douglasgscofield/bioinfo/tree/master/stacksTemplates)

Create templates for high-throughput genotyping platforms (iPlex, etc) from Stacks consensus sequences.


### [scripts](https://github.com/douglasgscofield/bioinfo/tree/master/scripts)

A handful of bioinformatics scripts:

**plink-pairwise-loci.pl**
 takes a [PLINK][] `.ped`-format input file and produce a matrix matching the `.mibs`-format pairwise output describing the number of complete loci (no missing data) contributing to each pairwise inbreeding coefficient in the `.mibs` file.

**fai2Bed.pl**
 creates BED intervals for each Fasta sequence described by the given Fasta index file.  Might be useful for separating input by chromosome for split-process-reduce parallelism.

**fai2WindowBed.pl**
 creates evenly-spaced BED intervals along Fasta sequences described by the given Fasta index file.  This is useful for creating BED files to guide sliding-window analyses with other tools.

**fillBed**
 fills BED-specified intervals of Fasta sequences with a given single-character sequence.  Especially useful for converting BED intervals to `N`.  Requires BioPerl.

**fastaSort**
 sorts a file of Fasta sequences by identifier name, naturally, so that `ìd_1` is followed by `id_2`, rather than `id_10`.  Requires BioPerl and the Perl module `Sort::Naturally`.

**gffSort**
 is a small Bash script that sorts a GFF file first by sequence name in column 1, and then by position numerically in column 4.  It assumes you might have generated this GFF with MAKER, so it removes the `###` lines between gene models if they are present.

**stacksExtractStats.pl**
reads a [Stacks][] log file and produce a table summarizing sample-specific statistics.  It currently prints number of RAD-Tags, number of stacks and mean stack coverage, and it is written so that other statistics can easily be harvested from the output.

**mummer2Vcf.pl**
reads file of SNPs and indels called by the [Mummer][] program `show-snps -T` and produce a VCF-ish file, collapsing consecutive indel characters into a single indel and adding the missing first base from indels by reading from the reference sequence file.  The format produced is close to compliant with [VCF][].  Requires BioPerl.

**subsampleReads.pl**
randomly selects a fraction of FastQ-format reads.

**phredDectector.pl**
attempts to determine the [Phred-scale quality encoding](http://en.wikipedia.org/wiki/FASTQ_format) used in a FastQ-format file.  If everything looks reasonable, it simply prints `33`, `64` or `59` (this last for older Solexa sequences) to `stdout`.  Does not require BioPerl.

**fermiExtractContigs.pl**
creates Fasta-format contigs from a [fermi][]-format `*.fq.gz` FastQ-like scaftig files.  Writes Fasta to `stdout`, giving each Fasta sequence its fermi sequence name and a description that includes sequence length, number of non-redundant reads that built the scaftig, and median coverage of non-redundant reads along the scafftig.  Requires BioPerl.

**fermiExtractContigs_simple.sh**
creates Fasta-format contigs from a [fermi][]-format `*.fq.gz` FastQ-like scaftig files.  Unlike **fermiExtractContigs.pl**, this does not use Perl.  Writes Fasta to `stdout`, giving each Fasta sequence its fermi sequence name and a description that includes sequence length and the number of non-redundant reads that built the scaftig.

**windowWig**
reads a data stream (for example, coverage values by position within reference sequences) and
produces a USCS [WIG][] file that summarizes median values within nonoverlapping windows. 

**intervalBed**
reads a data stream with reference-position optionally marked with boolean values (for example, 
presence-absence by position within reference sequences) and produces a [BED][] file 
describing intervals in which the values are true.

**samHeader2Bed.pl** 
reads a SAM header and produces [BED][] file(s) after applying a few filtering criteria.

**pos2bed**
reads a 2-column chromosome,position file -- for example, the first two columns of a VCF or GFF/GTF -- and
creates a BED file to stdout containing intervals that cover the chromosomes, with boundaries at
the positions.

**pileup2pro.pl**
reads `samtools mpileup` format and produces a profile file suitable for input to [mlRho][].

**mergePileupColumns**
merges columns from each BAM in multi-BAM `samtools mpileup` output into single columns.

**extractFastaSeqs.pl**
extracts named sequences from a FASTA file, or everything but.  Does not use BioPerl.

**extractFasta.pl**
extracts named FASTA sequences from a `makeblastdb -parse_seqids`-built blast database, and optionally provide a range for a subsequence.  Uses BioPerl.

**extractFasta**
extracts named FASTA sequences from a FASTA file indexed with BioPerl's `Bio::DB::Fasta`

**trimFastq.pl**
hard-trims a given amount from the 5' or 3' end (or both) from each read in a FastQ-format file, and optionally trims each read from the 3' end to a maximum length.

**shuffleFastq.pl** and **deshuffleFastq.pl**
convert FastQ-format files from separate read 1/read 2 files to interleaved and back, with some other handy options.

**compress_md5.sh**
compress standard input to a given file while simultaneously computing MD5 checksum for the uncompressed data; useful alone or with the shuffle routines above

**fastaOneline**
collapse a multiline Fasta sequence onto a single line; useful for `grep` etc.

**fastaGC.pl**
analyses GC content of Fasta-format sequences a few different ways.

**convertSequence.pl** and **convertAlignment.pl**
convert between any of the sequence and alignment formats supported by BioPerl, which must be available.  `convertAlignment.pl` can also convert aligned sequences to degapped unaligned sequences.  If you need to line-wrap a Fasta file, use

    convertSequence.pl -if fasta -of fasta < unwrapped.fa > wrapped.fa

**gmhmmp2Fasta.pl** and **gmhmmp2Table.pl**
extract Fasta sequences and a summary table from output produced by the ORF-finding tool [gmhmmp][].

**cutadaptReportScript.sh**
collects results from `*.cutReport` files produced by [cutadapt][] to quicklyl produce a table of adapter trimming results.


[PLINK]: https://www.cog-genomics.org/plink2
[WIG]:  http://genome.ucsc.edu/goldenPath/help/wiggle.html
[Stacks]:  http://creskolab.uoregon.edu/stacks/
[Mummer]:  http://mummer.sourceforge.net
[VCF]:  http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[fermi]:  https://github.com/lh3/fermi
[BED]:  http://genome.ucsc.edu/FAQ/FAQformat.html#format1
[mlRho]:  http://guanine.evolbio.mpg.de/mlRho
[velvet]: http://www.ebi.ac.uk/~zerbino/velvet
[gmhmmp]: http://www.genepro.com/Manuals/EuGM/EuGM_usage.aspx
[cutadapt]: https://code.google.com/p/cutadapt/

