#!/usr/bin/perl  -w

# This script takes two input files:1. A file with the first column containing names and allsubsequent columns containing data. 
# Example:names   col1   col2   col3...gen1    675    654    123...2. A file with name of columns you want to extract.
# This should not contain a header and each name of columnshould be in a separate line. 
# Example:
# col1
# col2
# coln

use warnings;
use strict;

print "USAGE: perl select_columns.pl <columnsData> <ids_file>\n";
print "Please check that your column names and ids are unique and identical\n";
print "example perl crop_columns.pl data_in_columns.txt column_names.txt\n";
print"_________________________ \n";

# dataset organized in columns with titles
chomp (my $infile = $ARGV[0]);
open (IN, $infile);

# a list of IDs, one per line
chomp (my $infile2 = $ARGV[1]);
open (IN2, $infile2);

my $outfile = $infile2 . ".selCol.tsv";
open (OUT, ">$outfile");

chomp(my $header = <IN>);
my @header = split(/\t/, $header);
my $num_items = @header;

# Assign an index to each column to be found
my @indices = ();

while (my $ID = <IN2>) {
	chomp($ID);
	my ($index) = grep {$header[$_] eq $ID } 0 .. $#header;
	push (@indices, $index);
}

print OUT $header[0];

foreach my $item(@indices){
	print OUT"\t", $header[$item];
}
print OUT"\n";

# Print columns at positions indicated in @indices
while (my $line = <IN>){
	chomp $line; 
	my @temp = split (/\t/, $line);

	print OUT $temp[0];

	foreach my $item(@indices) {
		print OUT "\t", $temp[$item];
	}
	print OUT "\n";
}

print "Your results are at: ", $outfile, "\n";

close IN;
close IN2;
close OUT;

exit;
