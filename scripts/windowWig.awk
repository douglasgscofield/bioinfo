#!/usr/bin/gawk -f
#
# Copyright (c) 2012 Douglas G. Scofield, Umeå Plant Sciences Centre, Umeå, Sweden
# douglas.scofield@plantphys.umu.se
# douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.
#
# Say you have a data stream (optionally containing a header line) with the 
# following format:
#
# reference1 1 value1
# reference1 2 value2
# ...
# reference1 1102 value2
# reference2 1 value1
# reference2 2 value2
# ...
#
# This script will cross each reference in column 1 in windowsize chunks based
# on the position in the 2nd column, computing the median of the values in
# column 3 within each window.  The median is calculated based on the number of
# values seen within each window, not in the size of the window; the median of
# a window containing a single value is that single value; windows which
# contain no values are reported to have a median of 0.
# 
# Caveats: Positions (column 2) must be sorted in increasing order within each
# reference (column 1).  Note that the positions within each reference are
# assumed to be monotonically increasing in steps of 1 starting from
# ref_start_pos, regardless of whether the data stream actually contains values
# for every position.  A value is printed for every window across a reference
# starting from ref_start_pos to the last reported position within the
# reference.  Every position, so defined, within every reference is guaranteed
# to be covered by a single reported window.
#
# I use some gawk extensions: 
#    delete cov;
#    asort()
#    maybe the PROCINFO["sorted_in"] mechanism??
# awk is so fast...
#
# CHANGELOG
# 2012-11-27 : generalize a bit
# 2012-11-26 : first version of the script
#
# TODO
# --- handle command-line arguments?  bet it's easy but awk and i, we are
#     just starting out...
# --- on option, skip 'comment lines' in input stream
# --- remove gawk dependencies?
# -x- generalize and put on github

# array sorting function: compare values for sort on numerical value in ascending order, 
# break ties by using index value
function cmp_numeric(i1, v1, i2, v2)
{
    return (v1 != v2) ? (v1 - v2) : (i1 - i2)  
}

# compute the median value of values within array; the rest of the arguments
# are for printing the verbose message
function compute_window_median(array, r, wbegin, wend, wsize, verb)
{
    n = asort(array);
    if (n < wsize) 
        if (verb) {
            print r ": window from " wbegin " to " wend " has less than " wsize " positions: " n > "/dev/stderr";
        }
    if (n == 1) {
        #print "index is " n > "/dev/stderr";
        ans = array[n];
    } else if ((n % 2) == 1) {
        #print "index is " int(n/2) > "/dev/stderr";
        ans = array[int(n / 2)];
    } else {
        #print "index is " (n/2) " and " (n/2 + 1) > "/dev/stderr";
        ans = (array[n / 2] + array[n / 2 + 1]) / 2;
    }
    return ans;
}


BEGIN {
    # parameters
    FS = "\t";    # input column separator
    ref_col = 1;  # column of reference, starting with 1 by awk convention
    pos_col = 2;  # column of position within reference (increasing within each reference)
    val_col = 3;  # column of value to windowize across positions within references
    header = 1;   # do we have header line(s) on the input to skip?
    ref_start_pos = 1; # the position at which references start, 1 by convention
    windowsize = 50;
    verbose = 0;

    # how we sort the values within each window
    PROCINFO["sorted_in"] = "cmp_numeric";

    # operational variables to track our position
    ref = "";     # current reference
    ibegin = -1;  # beginning of the current window
    iend = 0;     # the position immediately beyond the current window
    # didja notice? the current window is [ibegin, iend)
}

{
    if (header && NR <= header) {  # skip header line(s)
        next;
    }
    if (ibegin == -1) {  # initializing...
        ref = $(ref_col);
        n_pos = 0;
        ibegin = ref_start_pos;
        iend = ibegin + windowsize;
        print "fixedStep chrom=" ref " start=" ibegin " step=" windowsize " span=" windowsize;
    }

    if ($(ref_col) != ref || $(pos_col) >= iend) { # we're exiting the window, compute median, print it, start over

        medcov = compute_window_median(cov, ref, ibegin, iend, windowsize, verbose);
        print medcov;

        if ($(ref_col) != ref) {  # we're here because we hit a new reference
            if (verbose) {
                print ref ": moving to a new reference sequence " $(ref_col) > "/dev/stderr";
            }
            ref = $(ref_col);
            ibegin = ref_start_pos;
            print "fixedStep chrom=" ref " start=" ibegin " step=" windowsize " span=" windowsize;
        } else {
            ibegin = iend;
        }
        iend = ibegin + windowsize;

        # now check whether we should skip more windows...
        while ($(pos_col) >= iend) {
            if (verbose) {
                print ref ": skipping a 0-position window from " ibegin " to " iend > "/dev/stderr";
            }
            print 0;
            ibegin = iend;
            iend = ibegin + windowsize;
        }

        delete cov; # clear the array; split("", cov); is a portable alternative
        n_pos = 0;
    }

    # we're (now) within the window, add the median coverage value, go to the next data line
    ++n_pos;
    cov[n_pos] = $(val_col);
}

END {
    if (n_pos > 0) { # we're exiting the file and still need to compute median
        medcov = compute_window_median(cov, ref, ibegin, iend, windowsize, verbose);
        print medcov;
    }
    if (verbose) {
        print ref ": window from " ibegin " to " iend " was the last, final position " $(pos_col) ", end of input data" > "/dev/stderr";
    }
}
