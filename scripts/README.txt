Various small scripts.

samHeader2bed.pl  - Read a SAM header and produce a BED file containing reference
        sequences satisfying min/max length criteria.  Useful if you want to do 
        something like call SNPs in contigs >= 1kb:

        samtools view -H your.bam | samHeader2bed.pl -minlength 1000 - > min1000.bed
        samtools mpileup -l min1000.bed -u -f ref.fa your.bam | bcftools view ...


pileup2pro.pl  - Convert pileup to profile format as used for input to mlRho.  This
        is a format in which bases present at each position in a reference sequence
        are enumerated.  This script simply converts the format, so any filtering on 
        base/mapping quality, etc. that you may wish to do should be done when 
        generating the pileup:

        samtools mpileup -B -C50 -q1 -f ref.fa your.bam | pileup2pro.pl > mlRho-input.txt

        mlRho estimates population genetic parameters from NGS data sequenced from a 
        diploid genome.

        Haubold B, P Pfaffelhuber, and M Lynch. 2010. mlRho - a program for estimating
            the population mutation and recombination rates from shotgun-sequenced
            diploid genomes.  Molecular Ecology 19 Supplement 1:277-284.

            http://guanine.evolbio.mpg.de/mlRho

        It implements methods described in

        Lynch M. 2008. Estimation of nucleotide diversity, disequilibrium coefficients, 
            and mutation rates from high-coverage genome-sequencing projects. Molecular
            Biology and Evolution 25:2421-2431. 

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

