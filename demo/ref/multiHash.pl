#!/usr/bin/perl -w
use strict;
print "multi hash demo.\n";
my $multiHash={
	"jone" => {
		age => 47,
		eyes => "brown",
		weight => 186,
	},
	"mary" => {
		age => 23,
		eyes => "hazel",
		weight => 128,
	},#
	"bill" => {
		age => 35,
		eyes => "blue",
		weight => 157,	
	},
};#

print "bill's age:$multiHash->{'bill'}->{'age'}\n";
print "bill's eyes:$multiHash->{'bill'}->{'eyes'}\n";
print "bill's weight:$multiHash->{'bill'}->{'weight'}\n";

foreach my $people (keys $multiHash){
	print $people.":\n";
	print "$people 's age:$multiHash->{$people}->{'age'}\n";
	print "$people 's eyes:$multiHash->{$people}->{'eyes'}\n";
	print "$people 's weight:$multiHash->{$people}->{'weight'}\n";
	
}
