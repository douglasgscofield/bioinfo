#!/usr/bin/env perl

# Prepare StringTie-produced GTF for input as a GFF hints file to braker
# * Remove comment and blank lines
# * change 'StringTie' in 2nd column to 'b2h'
# * convert info field to GFF format: gene_id "MSTRG.1" to gene_id=MSTRG.1
# * add tag src=E to info

use strict;
use warnings;
use feature 'say';

sub info_to_gff($) {
    my $line = shift;
    my @line = split(';', $line);
    my @newline;
    foreach (@line) {
        s/^\s+//;
        s/\s+$//;
        #say "field='$_'";
        next if /^$/;
        my @field = split(/ /, $_, 2);
        $field[1] =~ s/^"(.*)"$/$1/;
        push @newline, "$field[0]=$field[1]";
    }
    #return "'".join(';', @newline)."'";
    return join(';', @newline);
}

while (<>) {
    chomp;
    next if /^\s*$/;  # blank line
    next if /^\s*#.*$/;  # comment line, with '#' as first nonblank character
    my @line = split /\t/;
    die "Line '$_' not produced by StringTie" if $line[1] ne "StringTie";
    $line[1] = "b2h";
    $line[8] = info_to_gff($line[8]) . ";src=E";
    say join("\t", @line);
}
