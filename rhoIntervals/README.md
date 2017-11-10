rhoIntervals : Infer recombination rates in arbitrary genomic intervals
============

Infer recombination rates (rho, &rho;) in arbitrary genomic intervals, using length-weighted contributions within underlying intervals.

Three files are expected for input, all with tab-separated columns:

* a four-column BED file describing recombination rate estimates in genomic intervals (`--rhofile`)
* a file specifying positions between which recombination rates should be inferred.  This file can have one of two formats, determined by the option used to specify the file
    * a two-column reference,position file describing genomic positions between which values of &rho; should be inferred (`--posfile`)  **OR**
    * a three-column BED file describing genomic intervals within which &rho; should be inferred (`--bedfile`)
* an FAI file describing the sizes of reference sequences named within the &rho; and position files (`--faifile`)


Options
-------

    --rhofile FILE    BED file of rho intervals
    --faifile FILE    File containing Fasta index (.fai) for reference

    --posfile FILE    File with SNP ref and position described in cols 1 and 2
    --bedfile FILE    BED file specifying genomic intervals

    --rho-column INT  Column of --rhofile containing the rho estimate, numbered from 1 [default 4]

    --help | -?       help


Example
-------

For the provided example input files, the following output is produced.

    $ ./rhoIntervals.pl --rho example_rho.bed --pos example_pos.pos --fai example.fai
    Chr21   1   443920  NA
    Chr21   443920  444023  0.00042253
    Chr21   444023  444035  0.00042253
    Chr21   444035  444088  0.00042253
    Chr21   444088  444127  0.00042253
    Chr21   444127  444129  0.00042253
    Chr21   444129  444139  0.00042253
    Chr21   444139  444192  0.00042253
    Chr21   444192  444196  0.00042253
    Chr21   444196  444276  0.0004149805625
    Chr21   444276  444282  0.000414885
    Chr21   444282  444313  0.000414885
    Chr21   444313  444321  0.000414885
    Chr21   444321  444337  0.000414885
    Chr21   444337  444338  0.000414885
    Chr21   444338  444350  0.000414885
    Chr21   444350  444381  0.000414885
    Chr21   444381  444409  0.000414885
    Chr21   444409  444414  0.000414885
    Chr21   444414  444423  0.000414885
    Chr21   444423  444533  0.000414885
    Chr21   444533  444554  0.000414885
    Chr21   444554  444621  0.000414885
    Chr21   444621  444647  0.000414885
    Chr21   444647  444661  0.000414885
    Chr21   444661  444667  0.000414885
    Chr21   444667  444714  0.000414885
    Chr21   444714  444734  0.000414885
    Chr21   444734  444744  0.000414885
    Chr21   444744  444765  0.000414885
    Chr21   444765  444776  0.000414885
    Chr21   444776  444780  0.00039214125
    Chr21   444780  444783  0.00038456
    Chr21   444783  444785  0.00038456
    Chr21   444785  444807  0.00038456
    Chr21   444807  444832  0.00038456
    Chr21   444832  444867  0.00038456
    Chr21   444867  444902  0.00038456
    Chr21   444902  444940  0.00038456
    Chr21   444940  445008  0.00038456
    Chr21   445008  445020  0.00038456
    Chr21   445020  445027  0.00038456
    Chr21   445027  445028  0.00038456
    Chr21   445028  445038  0.00038456
    Chr21   445038  445067  0.000281408965517241
    Chr21   445067  445068  0.000277725
    Chr21   445068  445097  0.000277725
    Chr21   445097  445103  0.000277725
    Chr21   445103  445114  0.000195211363636364
    Chr21   445114  445139  0.00018696
    Chr21   445139  445166  0.00018696
    Chr21   445166  445199  0.00018696
    Chr21   445199  445200  0.00018696
    Chr21   445200  445247  0.00018696
    Chr21   445247  445253  0.00018696
    Chr21   445253  445274  0.00018696
    Chr21   445274  445275  0.00018696
    Chr21   445275  445284  0.00018696
    Chr21   445284  445352  0.00018696
    Chr21   445352  445392  0.000102720975
    Chr21   445392  445400  0.000100561
    Chr21   445400  445437  0.000100561
    Chr21   445437  445447  0.000100561
    Chr21   445447  445454  0.000100561
    Chr21   445454  445458  0.000100561
    Chr21   445458  445504  0.000100561
    Chr21   445504  445509  0.000100561
    Chr21   445509  445574  0.000100561
    Chr21   445574  445657  0.000100561
    Chr21   445657  445674  0.000100561
    Chr21   445674  445688  0.000100561
    Chr21   445688  445694  0.000100561
    Chr21   445694  445713  0.000100561
    Chr21   445713  445732  0.000100561
    Chr21   445732  445757  0.000100561
    Chr21   445757  445811  0.000100561
    Chr21   445811  445838  0.000100561
    Chr21   445838  445847  0.000100561
    Chr21   445847  445854  0.000100561
    Chr21   445854  445881  0.000100561
    Chr21   445881  445902  0.000100561
    Chr21   445902  445969  0.000100561
    Chr21   445969  445971  0.000100561
    Chr21   445971  446012  0.000100561
    Chr21   446012  446069  0.000100561
