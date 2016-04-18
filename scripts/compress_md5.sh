#!/bin/bash

function usage() {
cat <<__USAGE__

USAGE:  $0  gz|bz2|none       filename
        $0  gz_self|bz2_self  filename

Saves stdin to filename (optionally compressed) while simultaneously
computing MD5 checksum on the uncompressed content, saved to filename.md5.

First argument is the compression format: gz bz2 none
Second argument is the output filename without any compression suffix

If the compression format is gz_self or bz2_self, use filename as input
instead of stdin, and remove filename if compression is successful.

Output is two files:

   filename     - compressed if requested, with suffix as appropriate
   filename.md5 - md5sum-format file with MD5 checksum for uncompressed
                  content of the file; filename present in filename.md5
                  is filename without any compression suffix

__USAGE__
exit 1
}

FORMAT="$1"
FILE="$2"
SLUG=$(echo "$FILE" | sed 's/\//\\\//g')
[[ "$FORMAT" && "$FILE" ]] || usage
SUFFIX=${FORMAT%_self}
[[ "$SUFFIX" = "none" ]] && SUFFIX="" || SUFFIX=".$SUFFIX"
MD5FILE="${FILE}.md5"
OUTFILE="${FILE}${SUFFIX}"


case $FORMAT in
gz)
    cat /dev/stdin | tee >(md5sum | sed "s/-/$SLUG/" > "$MD5FILE") | gzip -c > "$OUTFILE"
    ;;
gz_self)
    cat "$FILE" | tee >(md5sum | sed "s/-/$SLUG/" > "$MD5FILE") | gzip -c > "$OUTFILE" && rm -f "$FILE"
    ;;
bz2)
    cat /dev/stdin | tee >(md5sum | sed "s/-/$SLUG/" > "$MD5FILE") | bzip2 -c > "$OUTFILE"
    ;;
bz2_self)
    cat "$FILE" | tee >(md5sum | sed "s/-/$SLUG/" > "$MD5FILE") | bzip2 -c > "$OUTFILE" && rm -f "$FILE"
    ;;
none)
    cat /dev/stdin | tee >(md5sum | sed "s/-/$SLUG/" > "$MD5FILE") > "$OUTFILE"
    ;;
*)
    echo "Unknown compression format '$FORMAT', must be one of: gz bz2 none gz_self bz2_self"
    exit 1
    ;;
esac

