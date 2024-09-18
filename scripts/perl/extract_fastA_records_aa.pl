#!/usr/bin/perl -w

use strict;

use Bio::SeqIO;

my $usage    = "$0 <fastA_file> <outFile> <ids_file>\n";
my $infile   = shift or die $usage;
my $outfile  = shift or die $usage;
my $ids_list = shift or die $usage;
open (ID, $ids_list);
open (OUT, ">$outfile");



while (my $id = <ID>){
	chomp $id;
	
	my $seq_in  = Bio::SeqIO->new(
	                              -format   => 'fasta',
				      -file     => $infile,
				      -alphabet => 'protein'
				     );

	while (my $seq = $seq_in->next_seq) {

		my $target_id = $seq->id;
		if ($target_id =~ /$id/){

    			print OUT ">", $seq->id, "\n";	
    			print OUT $seq->seq(), "\n";
		}
	}
}

close ID;
close OUT;

exit;
                              
                            
