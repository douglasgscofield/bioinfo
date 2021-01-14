#!/bin/bash

# Process a collection of cutAdapt reports and produce a table showing
# the number of times each adapter was trimmed.  Table is produced
# to stdout.

# This script attempts to detect two different cutadapt output formats.  I'd
# been previously using cutadapt 1.3 for my pipeline, which produced the
# adapter's given name and cut results on a single line, version "pre":
#
#     Adapter 'Adapter_1_Illumina' (ACACTCTTTCCCTACACGACGCTGTTCCATCT), length 32, was trimmed 65 times.
#
# When I switched to cutadapt 1.8, I found these are now produced on 
# separate lines, which I am calling version "1.8":
#
#     === Adapter 'Adapter_1_Illumina' ===
#     
#     Sequence: ACACTCTTTCCCTACACGACGCTGTTCCATCT; Type: variable 5'/3'; Length: 32; Trimmed: 12 times.
#
# With paired-read unique dual UD adapters, these include 'First read: ' and 'Second read: ' in the headers,
# and do not include the single quotes around the adapter names.  I call this "UD" format:
#
#     === First read: Adapter Adapter_1_Illumina ===
#     === Second read: Adapter Adapter_1_Illumina ===
#     
#     Sequence: ACACTCTTTCCCTACACGACGCTGTTCCATCT; Type: variable 5'/3'; Length: 32; Trimmed: 12 times.
#

first=$1

all_reports=("$@")

TmpID=$$

# First, determine if it is cutadapt 1.8 or earlier

Version=

if grep -q '^Adapter.*, was trimmed ' "$first"
then
    Version="pre"
elif grep -q '^=== First read: ' "$first"
then
    Version="UD"
elif grep -q '^Sequence.*; Trimmed: ' "$first"
then
    Version="1.8"
else
    echo "Unknown cutadapt output format"
    exit 1
fi
echo "Reading cutadapt output format $Version"

function extract_pre()
{
    report=$1
    grep ', was trimmed ' $report | cut -d" " -f2,5,8 | sed -e "s/[',]//g" -e 's/ /\t/g'
}


function extract_18()
{
    report=$1
    while read line1; do
        read line2; read line3;
        echo "$line1 $line3" | cut -d' ' -f3,11,13 | sed -e "s/[';]//g" -e 's/ /\t/g'
    done < <(grep -A 2 '^=== Adapter ' $report | grep -v '^--$')
}

function extract_UD()
{
    report=$1
    while read line1; do
        read line2; read line3;
        echo "$line1 $line3" | cut -d' ' -f5,13,15 | sed -e "s/[';]//g" -e 's/ /\t/g'
    done < <(grep -P -A 2 '^=== (First|Second) read: Adapter ' $report | grep -v '^--$')
}

function extract()
{
    if [[ "$Version" = "pre" ]] ; then
        extract_pre $1
    elif [[ "$Version" = "1.8" ]] ; then
        extract_18 $1
    elif [[ "$Version" = "UD" ]] ; then
        extract_UD $1
    fi
}


for report in ${all_reports[@]} ; do

    SampleName=${report%.cutReport}

    SampleOut=$SampleName.$TmpID.tmp
    touch $SampleOut
    if [ "$report" = "$first" ] ; then
            echo -e "adapter\tlength\t$SampleName" > $SampleOut
            extract $report >> $SampleOut
            Tables="$SampleOut"
    else
            echo -e "$SampleName" > $SampleOut
            extract $report | cut -f3 >> $SampleOut
            Tables="$Tables $SampleOut"
    fi
done

paste $Tables && rm -f $Tables

