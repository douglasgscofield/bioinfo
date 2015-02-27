#!/bin/bash
# #!/bin/bash -x

TempID=$$
ScriptDir=..
Script=intervalBed

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


# -------------------------------- intervalBed tests


function runTest()
{
    ThisTest=$1
    TestOptions=$2
    Input=$3
    CompressedInput=$4
    echo -n "$ThisTest - $TestOptions ${Input}${CompressedInput}"
    if [[ $CompressedInput && ! -f $Input ]] ; then
        gzip -d -c < $CompressedInput > $Input
    fi
    Expected=${ThisTest}.expected
    Temp1=$ThisTest.$TempID.1
    Exec="$ScriptDir/$Script $TestOptions $Input"
    eval $Exec > $Temp1
    if diff $Temp1 $Expected ; then
        echo " - PASSED"
        RemoveTempFiles
        if [[ $CompressedInput ]] ; then
            rm -f $Input
        fi
    else
        echo " - FAILED"
    fi
    ClearVariables

}

runTest intervalBed_test_01 "col=4" intervalBed_test.input
runTest intervalBed_test_02 "col=4 strict=1" intervalBed_test.input
runTest intervalBed_test_03 "col=5" intervalBed_test.input
runTest intervalBed_test_04 "col=5 strict=1" intervalBed_test.input
runTest intervalBed_test_05 "crit_val=3" intervalBed_test.input
runTest intervalBed_test_06 "crit_val=19" intervalBed_test.input
runTest intervalBed_test_07 "col=0" intervalBed_test.input

