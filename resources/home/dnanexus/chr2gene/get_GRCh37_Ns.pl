#!/usr/bin/perl
## Name:   get_GRCh37_Ns.pl
## Author: Gang Wu
## Date:   March 9, 2011
## Purpose: get the location of sequence gaps "N" in GRCh37 reference genome
## Date source: ls /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/*.fa

use strict;
use warnings;

my $dir = "/nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/";

open(L,"$dir/fa.files") or die "cant open $dir/fa.files\n";
while (<L>) {
	chomp;
	my $f = $_;
	print "processing $f...\n";
	my $out = "$f.Ns";
	
	my ($lines,$prior) = (0,0);
	#Rationale: Y=[ATCG]
	# Case 1: YN, $prior=0,current=1    => print N, set $prior=1
	# Case 2: NN, $prior=1,current=1    => next unless it is the end of the whole seq.
	# Case 3: NY, $prior=1,current=0    => print N, set $prior=0
	open(F,$dir.$f) or die "cant open $dir$f\n";	
	open(OUT,">$out") or die "cant open $out to write\n";	
	while(<F>){
	  chomp;
	  if(!/^\>/) {
		for (my $i=0; $i<length($_); $i++) {
			my $base = substr($_,$i,1);
			my $pos = $lines * 60 + $i + 1;
			if (defined($base)) {
				if ($base eq "N" && $prior == 0) {
					#YN
					print OUT "$pos";
					$prior = 1;
				}elsif ($base eq "N" && $prior ==1){
					#NN
					print OUT ",$pos\n" if($i==(length($_)-1) && length($_) < 60); 
					#note the last line of seq in the fasta file is always less than 60
					#next base unless the end of seq
				}elsif ($base ne "N" && $prior ==1){
					#NY
					print OUT ",",$pos-1,"\n";
					$prior = 0;
				}else{
					#YY
					#next base
				}
			}
		}
		$lines++;
	  }
	}
	close F;
	close OUT;
}
close L;
