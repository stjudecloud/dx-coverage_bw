#!/usr/bin/perl
## Name:   get_GC.pl
## Author: Gang Wu
## Date:   Dec. 23, 2010
## Purpose: get the GC ratio of all genic regions in mm9 in order to correlate with NGS coverage
## Date source: ls *.fa|grep -v NT

use strict;
use warnings;


####+++++set the annotation of regions for each base++++++++
my @chrs = (1..22,"X","Y");
my $path = "../../chr2gene";
foreach my $chr (@chrs) {
	my $annot_file = "$path/chr$chr.gene.txt";
	my @base2region;
	print "processing chr$chr...\n";
	print "  set annotation for each base on chr$chr with $annot_file\n";
	open(EX, $annot_file) or die "cant open $annot_file...\n";
	while (<EX>) {
		chomp;
		my ($region,$rgType,$rgStart,$rgEnd) = split(/\t/);
		for (my $base=$rgStart; $base<=$rgEnd; $base++) {
			if (defined($base2region[$base])) { #when a base belongs to more than 2 regions in the same or different genes
				$base2region[$base] = $base2region[$base].",".$region;
			}else{
				$base2region[$base] = $region;
			}
		}			
	}
	close EX;

	##get each base from FASTA file
	my $dir = "./";
	my $f = "$chr.fa";
	print "processing $f...\n";
	
	my %counts;
	my $lines=0;
	open(F,$dir.$f) or die "cant open $dir$f\n";
	while(<F>){
	  chomp;
	  if(!/^\>/) {
		for (my $i=0; $i<length($_); $i++) {
			my $nuc = uc(substr($_,$i,1));
			my $base = $lines * 60 + $i + 1;
			if(defined($nuc) && defined($base2region[$base])){ 
				my @regions = split(/,/,$base2region[$base]);
				foreach my $region (@regions) {
					if ($nuc eq "G" ||  $nuc eq "C") {
						$counts{$region}->{'GC'}++;
						$counts{$region}->{'total'}++;
					#}elsif($nuc eq "A" ||  $nuc eq "T"){
					#	$counts{$region}->{'total'}++;
					}else{ # "N"
						$counts{$region}->{'total'}++;
					#	next;
					}
				}
			}

		}
		$lines++;
	  }
	}
	close F;
	my $out = "$f.all.GC";
	open(OUT,">$out") or die "cant open $out to write\n";	
	print OUT "#Gene|Accession|Strand|Region\t#GC\tTotal\tRatio\n";
	foreach my $region (sort keys %counts) {
		if (defined($counts{$region}->{'total'}) && $counts{$region}->{'total'}>0) {
			if (defined($counts{$region}->{'GC'})) {
				print OUT "$region\t" . $counts{$region}->{'GC'}.
							   "\t" . $counts{$region}->{'total'}.
							   "\t" . sprintf("%.2f",$counts{$region}->{'GC'}/$counts{$region}->{'total'});
				print OUT "\n";
			}else{
				print OUT "$region\t0\t" . $counts{$region}->{'total'}."\t0\n";
			}
		}
	}		
	close OUT;
}
