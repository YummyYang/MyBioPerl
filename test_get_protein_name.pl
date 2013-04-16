#! /usr/bin/perl -w
use strict;
use MyBioPerl;

info();
print "Print IUB table:\n";
print_iub_table();

print "\nOOPS, what I want is Amino Acid Table ...\n";
display_amino_acid();

print "\nand more : display amino acid, sort by 3.\n";
display_amino_acid_sort_by_3();

print "\nand more boring thing : display amino acid , sort by 1.\n";
display_amino_acid_sort_by_1();
