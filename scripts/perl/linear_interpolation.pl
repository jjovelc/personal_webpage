#!/usr/bin/perl -w

chomp (my $infile  = $ARGV[0]); # the file that comes out from MEGAN

chomp (my $interv  = $ARGV[1]); # a user-created input file
				# to create this file, open the
				# <infile> file in excel and find the 
				# library with the smallest number of
				# reads sampled. Copy the first number
				# of reads rarefied and estimate a number
				# by trial and error, that added to the first
				# number and iteratively to the resulting sub-total
				# nineteen times, will produce a number equal or 
				# smaller than the greatest number of rarefaction for
				# such library (en example follows)

=head1
Example

#Series: SFX3-0			
0	0		
5656	68.3	5656	5500
11312	69.1	11156	
16968	69.2	16656	
22624	69.2	22156	
28280	69.2	27656	
33936	69.3	33156	
39592	69.4	38656	
45248	69.4	44156	
50904	69.5	49656	
56560	69.6	55156	
62216	69.6	60656	
67872	69.6	66156	
73528	69.6	71656	
79184	69.6	77156	
84840	69.7	82656	
90496	69.7	88156	
96152	69.7	93656	
101808	69.8	99156	
107464	69.8	104656	
113120	69.8	110156	

=cut

my $outfile = $infile . ".xls";
open (IN1, $infile);
open (IN2, $interv);
open (OUT, ">$outfile");

# put intervals into an array
my @intervals;

while (my $int = <IN2>){
	push (@intervals, $int);
}
close IN2;	
	

my %megan_res;
do { 
	chomp (my $line = <IN1>);

	if ($line =~ m/#/){
		$line =~ s/#Series: //;
		print OUT "$line\n"; 
	} elsif ($line =~ m/^[0-9]/){
		

		$megan_res{'0'} = [0,0.0];

		for (my $i = 1; $i < 21; $i++){
			chomp(my $line1 = <IN1>);
			my @temp = split (/\t/, $line1);
			$megan_res{$i} = [$temp[0], $temp[1]];
		}

		foreach my $x2 (sort {$a <=> $b} @intervals){
			chomp $x2;
			foreach my $key (sort {$a <=> $b} keys %megan_res) {
				if ($megan_res{$key}[0] > $x2) {
					my $x3 = $megan_res{$key}[0];
					my $y3 = $megan_res{$key}[1];
					my $new_key = $key - 1;
					my $x1 = $megan_res{$new_key}[0];
					my $y1 = $megan_res{$new_key}[1];
					my $y2 = sprintf "%.2f", ((($x2 - $x1) * ($y3 - $y1)) / ($x3 - $x1)) + $y1;
					
					print OUT "$x2\t$y2\n"; 

					last;
				}
			}
		}


	}

} until eof;
		
close IN1;
close OUT;

exit;
