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

    $ ./rhoIntervals.pl --rhofile example_rho.bed --fai example.fai --pos example_pos.pos
    Chr21	0	443919	NA
    Chr21	443919	444022	0.00042253
    Chr21	444022	444034	0.00042253
    Chr21	444034	444087	0.00042253
    Chr21	444087	444126	0.00042253
    Chr21	444126	444128	0.00042253
    Chr21	444128	444138	0.00042253
    Chr21	444138	444191	0.00042253
    Chr21	444191	444195	0.00042253
    Chr21	444195	444275	0.000415076125
    Chr21	444275	444281	0.000414885
    Chr21	444281	444312	0.000414885
    Chr21	444312	444320	0.000414885
    Chr21	444320	444336	0.000414885
    Chr21	444336	444337	0.000414885
    Chr21	444337	444349	0.000414885
    Chr21	444349	444380	0.000414885
    Chr21	444380	444408	0.000414885
    Chr21	444408	444413	0.000414885
    Chr21	444413	444422	0.000414885
    Chr21	444422	444532	0.000414885
    Chr21	444532	444553	0.000414885
    Chr21	444553	444620	0.000414885
    Chr21	444620	444646	0.000414885
    Chr21	444646	444660	0.000414885
    Chr21	444660	444666	0.000414885
    Chr21	444666	444713	0.000414885
    Chr21	444713	444733	0.000414885
    Chr21	444733	444743	0.000414885
    Chr21	444743	444764	0.000414885
    Chr21	444764	444775	0.000414885
    Chr21	444775	444779	0.0003997225
    Chr21	444779	444782	0.00038456
    Chr21	444782	444784	0.00038456
    Chr21	444784	444806	0.00038456
    Chr21	444806	444831	0.00038456
    Chr21	444831	444866	0.00038456
    Chr21	444866	444901	0.00038456
    Chr21	444901	444939	0.00038456
    Chr21	444939	445007	0.00038456
    Chr21	445007	445019	0.00038456
    Chr21	445019	445026	0.00038456
    Chr21	445026	445027	0.00038456
    Chr21	445027	445037	0.00038456
    Chr21	445037	445066	0.000285092931034483
    Chr21	445066	445067	0.000277725
    Chr21	445067	445096	0.000277725
    Chr21	445096	445102	0.000277725
    Chr21	445102	445113	0.000203462727272727
    Chr21	445113	445138	0.00018696
    Chr21	445138	445165	0.00018696
    Chr21	445165	445198	0.00018696
    Chr21	445198	445199	0.00018696
    Chr21	445199	445246	0.00018696
    Chr21	445246	445252	0.00018696
    Chr21	445252	445273	0.00018696
    Chr21	445273	445274	0.00018696
    Chr21	445274	445283	0.00018696
    Chr21	445283	445351	0.00018696
    Chr21	445351	445391	0.00010488095
    Chr21	445391	445399	0.000100561
    Chr21	445399	445436	0.000100561
    Chr21	445436	445446	0.000100561
    Chr21	445446	445453	0.000100561
    Chr21	445453	445457	0.000100561
    Chr21	445457	445503	0.000100561
    Chr21	445503	445508	0.000100561
    Chr21	445508	445573	0.000100561
    Chr21	445573	445656	0.000100561
    Chr21	445656	445673	0.000100561
    Chr21	445673	445687	0.000100561
    Chr21	445687	445693	0.000100561
    Chr21	445693	445712	0.000100561
    Chr21	445712	445731	0.000100561
    Chr21	445731	445756	0.000100561
    Chr21	445756	445810	0.000100561
    Chr21	445810	445837	0.000100561
    Chr21	445837	445846	0.000100561
    Chr21	445846	445853	0.000100561
    Chr21	445853	445880	0.000100561
    Chr21	445880	445901	0.000100561
    Chr21	445901	445968	0.000100561
    Chr21	445968	445970	0.000100561
    Chr21	445970	446011	0.000100561
    Chr21	446011	446068	0.000100561
