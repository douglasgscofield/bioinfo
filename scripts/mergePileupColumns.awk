#!/usr/bin/awk -f

# Copyright (c) 2012 Douglas G. Scofield, Umeå Plant Sciences Centre, Umeå, Sweden
# douglas.scofield@plantphys.umu.se
# douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.
#
# Join mpileup columns into single columns when multiple BAMs were specified on the
# samtools mpileup command line.  Correctly handles 0-coverage BAMs.  See 
# https://github.com/douglasgscofield/bioinfo/tree/master/scripts#mergepileupcolumnsawk

BEGIN {
    FS = "\t";
    OFS = "\t";
    n_bams = 0;
}
{
    # $1 ref
    # $2 cov
    # $3 refbase
    # $4+ coverage
    cov=0;
    base_call="";
    base_qual="";
    map_qual="";
    i = 4;
    n_bams_now = 0;
    while (i <= NF) {  # walk through the columns
        if ($(i) > 0) {  # there is coverage here, 3 or 4 columns to the right
            cov += $(i);
            base_call = base_call $(i+1);
            base_qual = base_qual $(i+2);
            if (mpileup_s) {
                map_qual = map_qual $(i+3);
                i += 4;
            } else {
                i += 3;
            }
            ++n_bams_now;
        } else {  # no coverage, 3 columns to the right with no info
            i += 3;
            ++n_bams_now;
        }
    }
    if (! n_bams) {
        n_bams = n_bams_now;
    }
    if (! base_call) {  # we found no coverage at all
        base_call = "*";
        base_qual = "*";
        map_qual = "*";  # not written unless mpileup_s nonzero
    }

    if (mpileup_s) {
        print $1, $2, $3, cov, base_call, base_qual, map_qual;
    } else {
        print $1, $2, $3, cov, base_call, base_qual;
    }

    if (n_bams_now != n_bams) {
        print "Whoa, I was expecting " n_bams " BAM files, but on line " NR " you show me pileup for " n_bams_now ", seriously?!?" > "/dev/stderr"
        exit 1;
    }
}

