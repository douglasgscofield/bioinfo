Bioinformatics tools
====================

More details are provided within each subfolder.

[pileVar](https://github.com/douglasgscofield/bioinfo/tree/master/pileVar)
-------

`pileVar.pl` identifies variants based on mpileup output from
[samtools](http://samtools.sourcefourge.net).  This tool is currently geared
toward identifying fixed differences and regions of ambiguity due to indels.


[scripts](https://github.com/douglasgscofield/bioinfo/tree/master/scripts)
-------

A handful of bioinformatics scripts:

`windowWig`
reads a data stream (for example, coverage values by position within reference sequences) and
produces a USCS [WIG][] file that summarizes median values within nonoverlapping windows. 

`intervalBed`
reads a data stream with reference-position marked boolean values (for example, 
presence-absence by position within reference sequences) and produces a [BED][] file 
describing intervals in which the values are true.

`samHeader2Bed.pl` 
reads a SAM header and produces [BED][] file(s) after applying a few filtering criteria.

`pileup2pro.pl`
reads `samtools mpileup` format and produces a profile file suitable for input to [mlRho][].

`mergePileupColumns`
merges columns from each BAM in multi-BAM `samtools mpileup` output into single columns.

`extractFastaSeqs.pl`
extracts named sequences from a FASTA file, or everything but.

`shuffleFastq.pl` and `deshuffleFastq.pl`
convert FastQ-format files from separate read 1/read 2 files to interleaved and back.  These are based on similar scripts provided with [velvet][].


[WIG]:  http://genome.ucsc.edu/goldenPath/help/wiggle.html
[BED]:  http://genome.ucsc.edu/FAQ/FAQformat.html#format1
[mlRho]:  http://guanine.evolbio.mpg.de/mlRho
[velvet]: http://www.ebi.ac.uk/~zerbino/velvet
