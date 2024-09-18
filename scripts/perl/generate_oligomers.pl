#!/usr/bin/perl -w

use strict;
open (OUT, ">test.txt");

sub append_base {
    my ($length, $seq, $bases) = @_;
    
    if (length($seq) >= $length) {
    	print OUT "$seq\n";
    }
    else {
	append_base ($length, $seq . $_, $bases) for @$bases;
    }
}
				   
my @bases = ('A','C','G','T');
append_base (12, '', \@bases);
