#!/bin/bash
# #!/bin/bash -x

TempID=$$
ScriptDir=..
ScriptName=pos2bed

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


# -------------------------------- pos2bed_test_01


ThisTest="pos2bed_test_01"
echo -n $ThisTest - general test of script operation with a single chromosome
Options=${ThisTest}.options
Expected=${ThisTest}.expected
Temp1=$ThisTest.$TempID.1
Temp2=$ThisTest.$TempID.2
$ScriptDir/$ScriptName $(< $Options) > $Temp1
if diff $Temp1 $Expected ; then
    echo " - PASSED"
    RemoveTempFiles
else
    echo " - FAILED"
fi
ClearVariables

