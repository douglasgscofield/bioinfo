#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $o_minlen = 0;
my $o_single = "";
my $o_simplify = 0;
my $o_stdin = 0;
my $o_md5 = 0;
my $compress_md5 = "compress_md5.sh";
my ($filename1, $filename2, $filenameIn);

my $usage = "
$0  [ options ] [ interleaved-file | - ]  read1-outfile  read2-outfile

All files are FastQ format. Filenames that end with .gz or .bz2 are assumed to
be gzip- or bzip2-compressed, respectively.  Results with ill-formed FastQ files
are undefined.

  --minlen INT   Minimum length of sequence to keep.  Both members of a pair 
                 are dropped if one is too short [default: $o_minlen]
  --single FILE  Store reads which meet the --minlen criterion but have a mate
                 that did not meet it into this file
  --simplify     Output a single '+' for each third FastQ line in all outfiles,
                 removing any redundant (and wasteful) copy of the read name
  --md5          Calculate md5 checksums for the uncompressed FastQ stream for
                 each of the outfiles, saved to outfile.md5 after any .gz or .bz2
                 extension is removed.  This requires the $compress_md5 script
                 (available at https://github.com/douglasgscofield/bioinfo) to be
                 in your PATH.
  -              Read uncompressed FastQ from STDIN

";

GetOptions("minlen=i" => \$o_minlen, 
           "single=s" => \$o_single,
           "simplify" => \$o_simplify,
           "md5"      => \$o_md5,
           ""         => \$o_stdin) or die "$usage";

die "$usage" if (($#ARGV + 1) != (2 + !$o_stdin)) or $o_minlen < 0;
if ($o_single and !$o_minlen) {
    print STDERR "$0: no positive --minlen so no '$o_single' file created\n" if $o_single and !$o_minlen;
    $o_single = "";
}
if ($o_stdin) {
    $filenameIn = "/dev/stdin";
    $filename1  = $ARGV[0];
    $filename2  = $ARGV[1];
} else {
    $filenameIn = $ARGV[0];
    $filename1  = $ARGV[1];
    $filename2  = $ARGV[2];
}

sub md5_cmd($) {
    # create compress_md5.sh command
    my $filename = shift;
    my $format = "none";
    if      ($filename =~ s/\.gz$//) {
        $format = "gz";
    } elsif ($filename =~ s/\.bz2$//) {
        $format = "bz2";
    }
    my $cmd = "$compress_md5 $format '$filename'";
    return $cmd;
}

sub file_open {
    # open for reading or writing to/from pipe, with decompressing/compressing
    my ($filename, $mode, $do_md5) = @_;
    die "MD5 only available with mode '>'" if ($do_md5 && $mode ne ">");
    $do_md5 = md5_cmd($filename) if $do_md5;
    my ($comp_option, $cat_option);
    if ($mode eq "<") {
        $mode = "-|";  $comp_option = "-d"; $cat_option =""; # decompress file into pipe
    } elsif ($mode eq ">") {
        $mode = "|-";  $comp_option = ">"; $cat_option =">"; # compress from pipe into file
    } else {
        die "invalid mode $mode";
    }
    my ($fd, $fd_ok);
    if ($do_md5) {
    	$fd_ok = open($fd, $mode, "$do_md5");
    } elsif ($filename =~ /\.gz$/) {
    	$fd_ok = open($fd, $mode, "gzip -c $comp_option ${filename}");
    } elsif ($filename =~ /\.bz2$/) {
    	$fd_ok = open($fd, $mode, "bzip2 -c $comp_option ${filename}");
    } else {
    	$fd_ok = open($fd, $mode, "cat $cat_option ${filename}");
    }
    die "could not execute '$do_md5'" if ($do_md5 and not $fd_ok);
    die "could not open $filename" if not $fd_ok;
    return $fd;
}

my $INFILE   = file_open($filenameIn, "<");
my $OUTFILE1 = file_open($filename1,  ">", $o_md5);
my $OUTFILE2 = file_open($filename2,  ">", $o_md5);
my $SINGLE   = file_open($o_single,   ">", $o_md5) if $o_single;

while(<$INFILE>) {
    my $f1l1 = $_;
    my $f1l2 = <$INFILE>;
    my $f1l3 = <$INFILE>;
    my $f1l4 = <$INFILE>;

    my $f2l1 = <$INFILE>;
    my $f2l2 = <$INFILE>;
    my $f2l3 = <$INFILE>;
    my $f2l4 = <$INFILE>;

    $f1l3 = $f2l3 = "+\n" if $o_simplify;

    if ($o_minlen) {
        my $s1 = $f1l2;  # read 1 sequence
        chomp $s1;
        my $s2 = $f2l2;  # read 2 sequence
        chomp $s2;
        if (length($s1) >= $o_minlen && length($s2) >= $o_minlen) {
            print $OUTFILE1 $f1l1, $f1l2, $f1l3, $f1l4;
            print $OUTFILE2 $f2l1, $f2l2, $f2l3, $f2l4;
        } elsif ($o_single) {
            if (length($s1) >= $o_minlen) {
                print $SINGLE $f1l1, $f1l2, $f1l3, $f1l4;
            } elsif (length($s2) >= $o_minlen) {
                print $SINGLE $f2l1, $f2l2, $f2l3, $f2l4;
            }
        }
    } else {
        print $OUTFILE1 $f1l1, $f1l2, $f1l3, $f1l4;
        print $OUTFILE2 $f2l1, $f2l2, $f2l3, $f2l4;
    }
}

