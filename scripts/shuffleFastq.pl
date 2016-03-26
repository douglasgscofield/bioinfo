#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $o_stdout = 0;
my $o_subset = 0;
my $o_md5 = 0;
my $compress_md5 = "compress_md5.sh";
my $o_verbose = 0;
my $filename1;
my $filename2;
my $filenameOut;

my $usage = "
$0  [ options ]  read1-infile  read2-infile [ interleaved-outfile | - ]

All files are FastQ format. Filenames that end with .gz or .bz2 are assumed to
be gzip- or bzip2-compressed, respectively. Results with ill-formed FastQ files
are undefined.

  --subset INT   write INT read pairs, then stop
  --md5          Calculate md5 checksums for the uncompressed FastQ stream for
                 the outfile, saved to outfile.md5 after any .gz or .bz2
                 extension is removed.  This requires the $compress_md5 script
                 (available at https://github.com/douglasgscofield/bioinfo) to be
                 in your PATH.  This cannot be combined with '-'.
  -              write output to STDOUT

";

GetOptions("subset=i" => \$o_subset,
           "md5"      => \$o_md5,
           ""         => \$o_stdout) or die "$usage";
die "$usage" if (scalar(@ARGV) != (2 + !$o_stdout)) or ($o_md5 and $o_stdout);

$filename1   = $ARGV[0];
$filename2   = $ARGV[1];
$filenameOut = $o_stdout ? "/dev/stdout" : $ARGV[2];

# md5_cmd() and file_open() should be identical in shuffleFastq.pl and
# deshuffleFastq.pl

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

my $FILE1   = file_open($filename1,   "<");
my $FILE2   = file_open($filename2,   "<");
my $OUTFILE = file_open($filenameOut, ">", $o_md5);

while (<$FILE1>) {
	my $f1l1 = $_;
    my $f1l2 = <$FILE1>;
    my $f1l3 = <$FILE1>;
    my $f1l4 = <$FILE1>;
    print $OUTFILE $f1l1, $f1l2, $f1l3, $f1l4;
    my $f2l1 = <$FILE2>;
    my $f2l2 = <$FILE2>;
    my $f2l3 = <$FILE2>;
    my $f2l4 = <$FILE2>;
    print $OUTFILE $f2l1, $f2l2, $f2l3, $f2l4;

    if ($o_subset) {
		--$o_subset;
		print STDERR "subset $o_subset\n" if $o_verbose and ($o_subset % 100000 == 0);
		last if $o_subset <= 0;
	}
}

