Bioinformatics tools
====================

Each folder has its own README.md file with more details.

[pileVar](https://github.com/douglasgscofield/bioinfo/tree/master/pileVar)
-------

`pileVar.pl` identifies variants based on mpileup output from
[samtools](http://samtools.sourcefourge.net).  This tool is currently geared
toward identifying fixed differences and regions of ambiguity due to indels.


[scripts](https://github.com/douglasgscofield/bioinfo/tree/master/scripts)
-------

Various small bioinformatics scripts:

**`samHeader2bed.pl`** reads a SAM header and produces BED file(s) after applying a few filtering criteria.

**`pileup2pro.pl`** reads `samtools mpileup` format and produces a profile file suitable for input to
[mlRho][].

**`extractFastaSeqs.pl`** extracts named sequences from a FASTA file, or everything but.

**`shuffleFastq.pl`** and **`deshuffleFastq.pl`** convert FastQ-format files from separate read 1/read 2 files to interleaved and back.  These are based on similar scripts provided with [velvet][].

[mlRho]:  http://guanine.evolbio.mpg.de/mlRho

[velvet]: http://www.ebi.ac.uk/~zerbino/velvet
