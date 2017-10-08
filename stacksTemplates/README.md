Produce SNP-assay templates from Stacks batch files
===================================================

Read [Stacks](http://catchenlab.life.illinois.edu/stacks/) `batch_` files
(`.tags.tsv.gz`, `.snps.tsv.gz` and `.sumstats.tsv`) to produce templates for
high-throughput genotyping.  Within locus consensus sequences, template
candidates are filtered and central SNPs are identified based on minimum flank
sizes and other options.  Consensus sequences containing central SNPs that
fulfill a few quality criteria are then used to create sequences in which the
central SNP is called out and SNPs in flanking regions are assigned IUPAC
codes.

Example output:

~~~~
locus_0005448	TAACAGCTGGCTCCTCGGCGCTCAACCATCAACGTCCATCATCTCCGAGT[C/T]CTCCACCTCAGTATACTCCTCCCCGTCCTCAATCCTCGTACCCTCGGCAT
locus_0007269	AACATGTTGCTTTCCCCTCCTCTTTTATTTATGCAGCTTGACAAATGGCT[A/G]CACATAATGGTGTCAATTCGCCTCGATAATTTCACTTTGAGCATGCGGGC
locus_0007492	AATGAAGGAGAAAGGGGAGGAGGAAGACTTGTAGGTCACAGCCACAGCCA[A/G]GGTAAGTCCAAAAAATTGTTCTCTCGTTTTTTTATCTTTTTCRTTTTCGT
locus_0010785	GAAATGAATATGTCTTTGTTCTTCTCCATGGTTTGCTCATGAWTCATCTC[A/G]TTCTCACTCACTCACTCGCACACGCACACAAGRTGTTTGCTCYGTGTTAG
locus_0014827	GCAAGCTAAATAATTCTTTTAGGCAAACAAGTTTTGGTGTTCCTCCGCGC[A/G]GAGGGCCAAAGYGAGAGCGCGAGGCGAATAATTMTTGCAAATCACAAGAA
~~~~

A variety of options are designed to filter based on criteria for a particular
project in Sweden.  The code is documented sufficiently that it should be
possible to apply other custom criteria.



Options
-------

Input files are given with these options.  The --dir option can be used alone.
    
    --dir DIR              directory containing batch_1.{catalog.{snps,tags}.tsv.gz,sumstats.tsv} files
    --snps FILE            .snps.tsv.gz file
    --tags FILE            .tags.tsv.gz file
    --sumstats FILE        .sumstats.tsv file

To be considered, stack consensus sequences must pass these basic criteria:

    --minflank INT         minimum flank length; central SNPs are those between flanks [50]
    --maxflanksnps INT     maximum number of SNPs in each flank [3]
    --centralsnpgap INT    minimum gap (bp) between central SNPs [5]

Only central SNPs that pass simple filtering according to the following criteria are considered:

    --crit_n_region INT    min number of regions stack must be scored in [2]
    --crit_nwithq INT      min number of populations with the minor allele [1]
    --crit_qfreq FLOAT     min minor allele frequency across all scored populations [0.01]

For all SNP minor alleles, binomial confints are calculated as well:

    --alpha FLOAT          alpha for minor allele frequency confint (Wilson's method) [0.05]

All central SNPs that are considered are then ranked using these criteria, in order:

    Highest number of populations with minor allele across all regions
    Highest number of regions with minor allele
    Highest number of regions stack was scored in
    Highest lower bound of minor allele frequency confint across all regions
    Highest minor allele frequency across all regions

The central SNP that is ranked most highly is used, and the others are shifted to the flanks.

These options control output.

    --name STRING          name prefix for templates, only one of --name or --extended allowed
    --extended STRING      use extended template names, with STRING for the nonvariable prefix
                           Extended names are particular to this study.
                           STRINGR#cP##p###_locusnumber0padded
                           STRING: the prefix
                           R#c   : the regions containing the minor allele
                                   R3A : all regions
                                   R1N : only north
                                   R1T : only other
                                   R1L : only oland
                                   R2X : north and other
                                   R2Y : north and oland
                                   R2Z : other and oland
                           P##   : Total number of populations with minor allele, 0-padded
                           p###  : Number north, other, oland populations containing minor allele
    --padding INT          0-padded with for number suffix on template name [7]
    --trimtemplate         trim templates to have minflank+SNP+minflank bases [1]
    --annotatetemplate     add annotation columns to template output, for filtering [1]
    --freqdigits INT       digits right of decimal point for allele frequences [4]
    --D_limit INT          only read this number of templates; for debugging
    --verbose              verbose debugging output

