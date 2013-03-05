#!/usr/bin/perl -w
use strict;
print "dynamic array : a dynamic array demo.\n";

# 1. at first , we need a 'array size'. 
#    in practical, this number is getting from a dynamic variable.	
my $array_size = 6;

my @array = ( '0,' x ($array_size-1) . "0");
print "@array\n";
