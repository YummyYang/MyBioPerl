#! /usr/bin/perl -w
use strict;
print "qw instead.\n";

my @array = qw{aa,bb,cc};
print @array;
print "\n";

@array = ("aaa","bbb","ccc");
print @array;
print "\n";

# don't need use \/
@array = ("/home/yangming1/code/autorun","bbb","ccc");
for(my $i = 0 ; $i < @array ; $i++){
	print $array[$i]." ";
}
print "\n";
