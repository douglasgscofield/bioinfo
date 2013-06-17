#!/bin/bash

# Process a collection of cutAdapt reports and produce a table showing
# the number of times each adapter was trimmed.  Table is produced
# to stdout.

first=$1

all_reports=("$@")

TmpID=$$

for report in ${all_reports[@]} ; do

    SampleName=${report%.cutReport}

    SampleOut=$SampleName.$TmpID.tmp
    touch $SampleOut
    if [ "$report" = "$first" ] ; then
            echo -e "adapter\tlength\t$SampleName" > $SampleOut
            grep ", was trimmed " $report | cut -d" " -f2,5,8 | sed -e "s/[',]//g" -e 's/ /\t/g' >> $SampleOut
            Tables="$SampleOut"
    else
            echo -e "$SampleName" > $SampleOut
            grep ", was trimmed " $report | cut -d" " -f8 | sed -e "s/[',]//g" -e 's/ /\t/g' >> $SampleOut
            Tables="$Tables $SampleOut"
    fi
done

paste $Tables && rm -f $Tables

