#!/bin/bash

# all tests assume -s provided to samtools mpileup

ScriptDir=..
Script=mergePileupColumns

ScriptOptions="-v mpileup_s=1"


$ScriptDir/$Script $ScriptOptions mergePileupColumns_test_01.input > mergePileupColumns_test_01.output.1 2> mergePileupColumns_test_01.output.2
if diff mergePileupColumns_test_01.output.1 mergePileupColumns_test_01.expected.1 > mergePileupColumns_test_01.diff.1 ; then
    if diff mergePileupColumns_test_01.output.2 mergePileupColumns_test_01.expected.2 > mergePileupColumns_test_01.diff.2 ; then
	    echo "mergePileupColumns_test_01 - PASSED"
	    rm -f mergePileupColumns_test_01.diff.1 mergePileupColumns_test_01.output.1 mergePileupColumns_test_01.diff.2 mergePileupColumns_test_01.output.2
    else
	    echo "mergePileupColumns_test_01 - FAILED"
	    echo "diff output in mergePileupColumns_test_01.diff.2"
    fi
else
    echo "mergePileupColumns_test_01 - FAILED"
    echo "diff output in mergePileupColumns_test_01.diff.1"
fi

$ScriptDir/$Script $ScriptOptions mergePileupColumns_test_02.input > mergePileupColumns_test_02.output.1 2> mergePileupColumns_test_02.output.2
if diff mergePileupColumns_test_02.output.1 mergePileupColumns_test_02.expected.1 > mergePileupColumns_test_02.diff.1 ; then
    if diff mergePileupColumns_test_02.output.2 mergePileupColumns_test_02.expected.2 > mergePileupColumns_test_02.diff.2 ; then
	    echo "mergePileupColumns_test_02 - PASSED"
	    rm -f mergePileupColumns_test_02.diff.1 mergePileupColumns_test_02.output.1 mergePileupColumns_test_02.diff.2 mergePileupColumns_test_02.output.2
    else
	    echo "mergePileupColumns_test_02 - FAILED"
	    echo "diff output in mergePileupColumns_test_02.diff.2"
    fi
else
    echo "mergePileupColumns_test_02 - FAILED"
    echo "diff output in mergePileupColumns_test_02.diff.1"
fi

$ScriptDir/$Script $ScriptOptions mergePileupColumns_test_03.input > mergePileupColumns_test_03.output.1 2> mergePileupColumns_test_03.output.2
if diff mergePileupColumns_test_03.output.1 mergePileupColumns_test_03.expected.1 > mergePileupColumns_test_03.diff.1 ; then
    if diff mergePileupColumns_test_03.output.2 mergePileupColumns_test_03.expected.2 > mergePileupColumns_test_03.diff.2 ; then
	    echo "mergePileupColumns_test_03 - PASSED"
	    rm -f mergePileupColumns_test_03.diff.1 mergePileupColumns_test_03.output.1 mergePileupColumns_test_03.diff.2 mergePileupColumns_test_03.output.2
    else
	    echo "mergePileupColumns_test_03 - FAILED"
	    echo "diff output in mergePileupColumns_test_03.diff.2"
    fi
else
    echo "mergePileupColumns_test_03 - FAILED"
    echo "diff output in mergePileupColumns_test_03.diff.1"
fi

$ScriptDir/$Script $ScriptOptions mergePileupColumns_test_04.input > mergePileupColumns_test_04.output.1 2> mergePileupColumns_test_04.output.2
if diff mergePileupColumns_test_04.output.1 mergePileupColumns_test_04.expected.1 > mergePileupColumns_test_04.diff.1 ; then
    if diff mergePileupColumns_test_04.output.2 mergePileupColumns_test_04.expected.2 > mergePileupColumns_test_04.diff.2 ; then
	    echo "mergePileupColumns_test_04 - PASSED"
	    rm -f mergePileupColumns_test_04.diff.1 mergePileupColumns_test_04.output.1 mergePileupColumns_test_04.diff.2 mergePileupColumns_test_04.output.2
    else
	    echo "mergePileupColumns_test_04 - FAILED"
	    echo "diff output in mergePileupColumns_test_04.diff.2"
    fi
else
    echo "mergePileupColumns_test_04 - FAILED"
    echo "diff output in mergePileupColumns_test_04.diff.1"
fi
