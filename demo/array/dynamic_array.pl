#!/usr/bin/perl -w
use strict;
print "dynamic array : a dynamic array demo.\n";

# 1. at first , we need a 'array size'. 
#    in practical, this number is getting from a dynamic variable.	
my $array_size = 6;

#my @array = ( '0,' x ($array_size-1) . "0"); 	# This is WRONG !!
#print "@array\n";

# this is not like the SQL codes below :
#my $query = "insert into $table ("
#	. join("," , @fieldnames)
#	. ") values ("
#	. "?, " x (@fieldnames - 1)  	# <-- not like this
#	. "?)";

# so we must initialize them in the normal way:
my @array;
$#array = $array_size -1;
foreach my $item (@array){
	$item = 1;
}

foreach my $item (@array){
	print "$item\n";
}

# or the C style normal way :

for(my $i=0; $i < $array_size ; $i++){
	$array[$i] = 2;
}

for(my $i=0 ; $i < $array_size ; $i++){
	print "array[$i] : $array[$i] \n";
}
