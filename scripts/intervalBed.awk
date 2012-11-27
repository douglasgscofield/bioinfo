#!/usr/bin/awk -f
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
#     reference1  1   0
#     reference1  2   1
#     reference1  3   1
#     reference1  4   1
#     reference1  5   1
#     reference1  6   1
#     reference1  7   1
#     reference1  8   1
#     reference1  9   0
#     reference1  10  0
#     reference1  11  1
#     reference1  12  1
#     reference1  13  1
#     reference1  14  0
#     reference1  15  0
#     reference1  16  0
#     ...
#
# This script will encode each stretch of 1s in the 3rd column into a BED 
# interval, for example the above input would result in the following BED
# lines:
#
#     reference1  1   8
#     reference1  10  13
#
# Note that intervals in BED are 0-based and [begin, end).  The script
# assumes that input positions are 1-based.  The script also assumes that
# positions not represented in the input have values of 0.
#
# Optionally, a grace distance may be defined such that stretches that are
# separated by less than or equal to the grace distance will be annotated
# as if they are connected.  For the above dataset, a grace greater than or
# equal to 2 will result in a single annotation line:
#
#     reference1  1   13
# 
# while a grace of 1 will not connect the stretches.  Grace distances
# also may apply across positions not represented in the input.  Grace
# only extends to intervals within the current reference.  Grace will only
# connect identifiable intervals and will never implicitly include the
# beginning or end of references.
#
# CHANGELOG
# 2012-11-27 : create the script
#
# TODO
# --- implement an inverse, connecting stretches of 0s rather than 1s,
#     where grace still applies and undefined positions on the input are 
#     still assumed to be 0
# --- handle command-line arguments?  bet it's easy but awk and i, we are
#     just starting out...
# --- on option, skip 'comment lines' in input stream
# --- skipping missing positions on input
# -x- monotonically increasing positions

BEGIN {
    # parameters
    FS = "\t";    # input column separator
    OFS = "\t";   # output column separator
    ref_col = 1;  # column of reference, starting with 1 by awk convention
    pos_col = 2;  # column of position within reference (increasing within each reference)
    val_col = 3;  # column of boolean values across positions within references
    header = 1;   # do we have header line(s) on the input to skip?
    ref_start_pos = 1; # the position at which references start, 1 by convention
    out_start_pos = 0; # the position at which references start, 1 by convention
    out_adjust = out_start_pos - ref_start_pos;  # add this to get output position
    grace = 50;   # grace distance, we also use a lookahead[] array
    inverse = 0;  # connect 0-valued intervals instead of 1-valued intervals; UNIMPLEMENTED
    verbose = 0;

    print_track = 1;
    trackname = "booleanIntervals";
    trackdesc = "intervals of " (inverse ? 0 : 1) "s, grace distance " grace;

    # operational variables to track our position
    ref = "";     # current reference
    prev_pos = -1;     # current position
    prev_val = -1;     # current value
    ibegin = -1;  # beginning of the current interval
    iend = 0;     # the position immediately beyond the current interval
    # didja notice? the current interval is [ibegin, iend)
    in_grace = 0; # counts the positions of grace we've been given 
    crit_val = 1; # what val must be for interval determination
    missing_val = 0; # what val is given to missing positions?
}

function print_interval(r, ib, ie)
{
    print r, (ib + out_adjust), (ie + out_adjust);
}

{
    if (header && NR <= header) {  # skip header line(s)
        next;
    }
    if (ibegin == -1) {  # initializing...
        ref = $(ref_col);
        pos = $(pos_col);
        val = $(val_col);
        ibegin = (val == crit_val) ? pos : 0;
        iend = ibegin + 1;
        if (print_track)
            print "track name=" trackname " description=\"" trackdesc "\"";

        next;
    }

    if ($(ref_col) != ref) {  # new reference seen
        if (ibegin) {  # we were working on an interval, ignore any grace
            print_interval(ref, ibegin, iend);
            ibegin = 0;
        }
        if (grace)
            in_grace = 0;
        ref = $(ref_col);
    }

    pos = $(pos_col);
    val = $(val_col);

    if (ibegin) {  # working on an interval
        if (val == crit_val) {
            if (pos == iend) { # we did not skip positions
                if (grace && in_grace) {  # use the grace we've been given
                    iend += in_grace;
                    in_grace = 0;
                }
                ++iend;
            } else if (pos > iend && grace && pos - iend <= grace) { 
                # we skipped positions on input, they were within the grace we have
                iend = pos + 1;
                in_grace = 0;
            } else {  # we have no grace, print previous interval and start a new one
                print_interval(ref, ibegin, iend);
                ibegin = pos;
                iend = ibegin + 1;
            }
        } else if (grace && ++in_grace <= grace) {
            # keep moving until no more grace
        } else {  # not crit_val and no grace left, print interval and reset
            if (grace)  # we cannot use grace here
                in_grace = 0;
            print_interval(ref, ibegin, iend);
            ibegin = 0; 
        }
    } else if (val == crit_val) {  # begin an interval
        ibegin = pos;
        iend = ibegin + 1;
    }
}

END {
    if (ibegin) { # working on an interval, ignore any grace
        print_interval(ref, ibegin, iend);
    }
}
