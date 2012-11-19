#! /usr/bin/perl -w
use strict;
print "pass 2 values to a function.\n";

#hash
my %hash = (
	"key1" => "vaule1",
	"key2" => "vaule2",
);
# ref to hash
my $ref=\%hash;

# a ref to anonymous hash
my %hash_ref = {
	"key1" => "vaule1",
	"key2" => "vaule2",
};

foreach (reverse keys %hash){
	print "$_ => $hash{$_}\n";
}

print "key1:".$hash{key1}."\n";
print "key2:".%$ref->{key2}."\n";
print "key2:".$ref->{key2}."\n";	# same.
