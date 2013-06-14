#!/usr/bin/perl

# Process gmhmmp output into a table of ORF discoveries

print "contig\torfnum\tstrand\tstart\tend\tlength\tclass\n";

use strict;
use warnings;

my $thiscontig;

while (<>) {
	next if ! /^(  *[0-9][0-9]*|FASTA)/;
	chomp;
	s/ +/ /;
	s/^ //;
	if (/^FASTA/) { # begins with "FASTA definition line: "
		$thiscontig = substr($_, 23);
	} else {
		my @orf = split ' ';
		unshift @orf, $thiscontig;
		print join ("\t", @orf), "\n";
	}
}

