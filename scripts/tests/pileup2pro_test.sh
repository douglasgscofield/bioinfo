#!/bin/bash
# #!/bin/bash -x

TempID=$$
ScriptDir=..

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


# -------------------------------- pileup2pro_test_01


ThisTest="pileup2pro_test_01"
echo $ThisTest - general test of script operation 
Input=pileup2pro_test_01.input
if [ "$CompressedInput" != "" -a ! -f $Input ] ; then
    gzip -d -c < $CompressedInput > $Input
fi
Options=pileup2pro_test_01.options
Expected=pileup2pro_test_01.expected
Temp1=$ThisTest.$TempID.1
Temp2=$ThisTest.$TempID.2
$ScriptDir/pileup2pro.pl $(< $Options) > $Temp1
if diff $Temp1 $Expected ; then
    echo "$ThisTest - PASSED"
    RemoveTempFiles
    if [ "$CompressedInput" != "" ] ; then
        rm -f $Input
    fi
else
    echo "$ThisTest - FAILED"
fi
ClearVariables


# -------------------------------- pileup2pro_test_02


ThisTest="pileup2pro_test_02"
echo $ThisTest - general test of script operation with restrictive --which-bams option
Input=pileup2pro_test_02.input
if [ "$CompressedInput" != "" -a ! -f $Input ] ; then
    gzip -d -c < $CompressedInput > $Input
fi
Options=pileup2pro_test_02.options
Expected=pileup2pro_test_02.expected
Temp1=$ThisTest.$TempID.1
Temp2=$ThisTest.$TempID.2
$ScriptDir/pileup2pro.pl $(< $Options) > $Temp1
if diff $Temp1 $Expected ; then
    echo "$ThisTest - PASSED"
    RemoveTempFiles
    if [ "$CompressedInput" != "" ] ; then
        rm -f $Input
    fi
else
    echo "$ThisTest - FAILED"
fi
ClearVariables


# -------------------------------- pileup2pro_test_03


ThisTest="pileup2pro_test_03"
echo $ThisTest - general test of script operation with --has-mapping-quality option
Input=pileup2pro_test_03.input
if [ "$CompressedInput" != "" -a ! -f $Input ] ; then
    gzip -d -c < $CompressedInput > $Input
fi
Options=pileup2pro_test_03.options
Expected=pileup2pro_test_03.expected
Temp1=$ThisTest.$TempID.1
Temp2=$ThisTest.$TempID.2
$ScriptDir/pileup2pro.pl $(< $Options) > $Temp1
if diff $Temp1 $Expected ; then
    echo "$ThisTest - PASSED"
    RemoveTempFiles
    if [ "$CompressedInput" != "" ] ; then
        rm -f $Input
    fi
else
    echo "$ThisTest - FAILED"
fi
ClearVariables


# -------------------------------- pileup2pro_test_04


ThisTest="pileup2pro_test_04"
echo $ThisTest - general test of script operation with --has-mapping-quality and --which-bams options
Input=pileup2pro_test_04.input
if [ "$CompressedInput" != "" -a ! -f $Input ] ; then
    gzip -d -c < $CompressedInput > $Input
fi
Options=pileup2pro_test_04.options
Expected=pileup2pro_test_04.expected
Temp1=$ThisTest.$TempID.1
Temp2=$ThisTest.$TempID.2
$ScriptDir/pileup2pro.pl $(< $Options) > $Temp1
if diff $Temp1 $Expected ; then
    echo "$ThisTest - PASSED"
    RemoveTempFiles
    if [ "$CompressedInput" != "" ] ; then
        rm -f $Input
    fi
else
    echo "$ThisTest - FAILED"
fi
ClearVariables


# -------------------------------- pileup2pro_test_05


ThisTest="pileup2pro_test_05"
echo $ThisTest - general test of script operation with --has-mapping-quality and restrictive --which-bams options
Input=pileup2pro_test_05.input
if [ "$CompressedInput" != "" -a ! -f $Input ] ; then
    gzip -d -c < $CompressedInput > $Input
fi
Options=pileup2pro_test_05.options
Expected=pileup2pro_test_05.expected
Temp1=$ThisTest.$TempID.1
Temp2=$ThisTest.$TempID.2
$ScriptDir/pileup2pro.pl $(< $Options) > $Temp1
if diff $Temp1 $Expected ; then
    echo "$ThisTest - PASSED"
    RemoveTempFiles
    if [ "$CompressedInput" != "" ] ; then
        rm -f $Input
    fi
else
    echo "$ThisTest - FAILED"
fi
ClearVariables

