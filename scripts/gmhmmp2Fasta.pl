#!/usr/bin/perl

# Process gmhmmp output into a Fasta file of ORFs

use strict;
use warnings;

my $infasta = 0;

while (<>) {
	chomp;
	$infasta = 1 if /^>gene/;
	if (/^$/) {
		print "\n" if $infasta == 1;
		$infasta = 0;
	}
	print "$_\n" if $infasta;
}
