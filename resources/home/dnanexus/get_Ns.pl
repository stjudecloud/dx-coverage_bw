#!/usr/bin/perl
## Name:   get_Ns.pl
## Author: Gang Wu
## Date:   March 9, 2011
## Purpose: get the location of sequence gaps "N" in GRCh37 reference genome
## Date source: ls /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/*.fa

use strict;
use warnings;

if  (@ARGV<1) {
   print "Usage: get_Ns.pl <path to chromosomes subdirectory of FASTA> <canonical_chr_sizes.txt>\n";
   print "\tWrites output to the current working directory\n"; 
   exit;
}

my $dir = $ARGV[0];
my $chrs = $ARGV[1]; 

# Open the list of canonical chromosomes
open(L,"$chrs") or die "cant open $chrs\n";
while (<L>) {
	chomp;
	my ($chr, $length) = split(/\t/, $_); 
	my $f = "$chr.fa";
	print "processing $f...\n";
	my $out = "$f.Ns";
	
	my ($lines,$prior, $pos) = (0,0,0);
	#Rationale: Y=[ATCG]
	# Case 1: YN, $prior=0,current=1    => print N, set $prior=1
	# Case 2: NN, $prior=1,current=1    => next unless it is the end of the whole seq.
	# Case 3: NY, $prior=1,current=0    => print N, set $prior=0
	open(F,$dir."/".$f) or die "cant open $dir$f\n";	
	open(OUT,">$out") or die "cant open $out to write\n";	
	while(<F>){
	  chomp;
	  if(!/^\>/) {
		my $length = length($_); 
		for (my $i=0; $i<length($_); $i++) {
			my $base = substr($_,$i,1);
			my $mpos = $pos + $i + 1;
			if (defined($base)) {
				if ($base eq "N" && $prior == 0) {
					#YN
					print OUT "$mpos";
					$prior = 1;
				}elsif ($base eq "N" && $prior ==1){
					#NN
					print OUT ",$mpos\n" if($i==(length($_)-1) && length($_) < 60); 
					#note the last line of seq in the fasta file is always less than 60 
					#next base unless the end of seq
				}elsif ($base ne "N" && $prior ==1){
					#NY
					print OUT ",",$mpos-1,"\n";
					$prior = 0;
				}else{
					#YY
					#next base
				}
			}
		}
		$pos += $length;
		$lines++;
	  }
	}
	close F;
	close OUT;
}
close L;
