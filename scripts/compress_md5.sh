#!/bin/bash

function usage() {
cat <<__USAGE__

USAGE:  $0  gz|bz2|none  filename

Saves stdin to filename (optionally compressed) while simultaneously
computing MD5 checksum on the uncompressed content, saved to filename.md5.

First argument is the compression format: gz bz2 none
Second argument is the output filename without any compression suffix

Output is two files:

   filename     - compressed if requested, with suffix as appropriate
   filename.md5 - md5sum-format file with MD5 checksum for uncompressed
                  content of the file; filename present in filename.md5
                  is filename without any compression suffix

__USAGE__
exit 1
}

FORMAT=$1
FILE=$2
[[ "$FORMAT" && "$FILE" ]] || usage
SUFFIX=$FORMAT
[[ "$SUFFIX" = "none" ]] && SUFFIX="" || SUFFIX=".$SUFFIX"
MD5FILE=${FILE}.md5
OUTFILE=${FILE}${SUFFIX}


case $FORMAT in
gz)
    cat /dev/stdin | tee >(md5sum | sed "s/-/$FILE/" > "$MD5FILE") | gzip -c > "$OUTFILE"
    ;;
bz2)
    cat /dev/stdin | tee >(md5sum | sed "s/-/$FILE/" > "$MD5FILE") | bzip2 -c > "$OUTFILE"
    ;;
none)
    cat /dev/stdin | tee >(md5sum | sed "s/-/$FILE/" > "$MD5FILE") > "$OUTFILE"
    ;;
*)
    echo "Unknown compression format '$FORMAT', must be one of: gz bz2 none"
    exit 1
    ;;
esac

