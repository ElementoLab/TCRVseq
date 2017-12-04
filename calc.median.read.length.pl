# !/usr/local/bin/perl -w
# Set use

use strict;
use warnings;

my $usage = "Usage: perl $0 <INFILE>";
my $infile = shift or die $usage;
open(IN, "<$infile") || die "Unable to open $infile: $!";
# Set variable for list of length values
my @lengths = ();
my $line_length = 0;
my $result = 0;
my $flag = 3;

while (my $line = <IN>) {
# Skip the sequence identifier lines
	if ($flag%4) {
	# Add  line lengths to array
		my $flag = $flag++;
		next;
		} else {
	# Add  line lengths to array
		chomp $line;
		my $line_length = length $line;
		push @lengths, $line_length;
		my $flag = $flag++;
		}
	}

my @sorted_lengths = sort {$a <=> $b} @lengths;

my $array_length = @sorted_lengths;
if ($array_length%2) {
	# Odd
	my $result = $sorted_lengths[($array_length/2)];
	print $result;
} else {
	# Even
	my $one = $sorted_lengths[($array_length/2)-1];
	my $two = $sorted_lengths[($array_length/2)];
	my $result = ($one + $two) / 2;
	print $result;
}

close(IN);
