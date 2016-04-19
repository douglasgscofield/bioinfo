#!/bin/bash

# extract contigs from fermi results file
#
# Fermi has a FastQ-like format, with scaftig name, scaftig length, and ?? number of non-duplicated reads?
#
# @72248153:1572501_0	136	23
# TGATTCAAGGCTATCTCCCCAGTGGAGTCGCCAGTTTCTGTTGTAGGAAAATCTGCCTGCGGACAAGTAGATCGACCTCTAAGGCTTGACCTTAACCCTAGACAGAAAAGTGCTTATCATGGACACTGAATTTCTT
# +
# "##$$%&&&&'(()*+,-./01122334445666788888888888888888888888888888888888888888888888888888888888888888876655433332110/.-,+*)((''&&%%%$###"

Scaftigs=$1

# grep out the scaftig names and the following line, which is the sequence
# convert the scaftig name to fasta-like, and the length and reads to part of the description

gzip -dc -f "$Scaftigs" | grep -A 1 '^@' | grep -v '^--$' | sed -e 's/^@/>/g' -e 's/\t\([0-9]*\)\t\([0-9]*\)/ length:\1,reads:\2/g'


