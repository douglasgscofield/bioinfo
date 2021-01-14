#!/bin/bash
# #!/bin/bash -x

TempID=$$
ScriptDir=..
ScriptName=braker_RNAseq_gtf.pl

# Complete set of variables used in tests; reset between tests

function ClearVariables()
{
    ThisTest=
    Input=
    CompressedInput=
    Expected=
    Output=
    Temp1=
    Temp2=
    Temp3=
}
ClearVariables

function RemoveTempFiles()
{
    rm -f $Temp1 $Temp2
    if [ "$Temp3" != "" ] ; then
        rm -f $Temp3
    fi
}
RemoveTempFiles


# -------------------------------- braker_RNAseq_gtf_test_01


ThisTest="braker_RNAseq_gtf_test_01"
echo -n $ThisTest - general test of script operation
Expected=${ThisTest}.expected
Input=${ThisTest}.input
Temp1=$ThisTest.$TempID.1
$ScriptDir/$ScriptName $Input > $Temp1
if diff $Temp1 $Expected ; then
    echo " - PASSED"
    RemoveTempFiles
else
    echo " - FAILED"
fi
ClearVariables

