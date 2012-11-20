#!/usr/bin/perl -w
use strict;
print "parse string to array.\n";

my $str="a:bb:ccc:dddd";
my @array = split /:/,$str;
print "array :@array\n";
