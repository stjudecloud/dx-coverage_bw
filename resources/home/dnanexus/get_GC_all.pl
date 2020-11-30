#!/usr/bin/perl
## Name:   get_GC_all.pl
## Author: Gang Wu
## Date:   Dec. 23, 2010
## Purpose: get the GC ratio of all genic regions in a genome in order to correlate with NGS coverage

use strict;
use warnings;
use Bio::Perl;
use Bio::SeqIO;

use File::Basename; 

if  (@ARGV<1) {
   print "Usage: get_GC_all.pl <path to output of chr2gene.pl> <FASTA file>\n";
   exit;
}


my ($path, $fasta_file) = @ARGV;

####+++++set the annotation of regions for each base++++++++
my $in = Bio::SeqIO->new (-file => $fasta_file,
                              -format => 'fasta');
while ( my $seq = $in->next_seq() ) {
	my $id = $seq->id;
	my $chr = $id;  
	$chr =~ s/^chr//;  
	my $annot_file = "$path/chr$chr.gene.txt";  
	next if (! -e $annot_file);
	my %counts;
	my $lines=0;
	my $chr_seq = $seq->seq;

	print "processing chr$chr...\n";
	open(EX, $annot_file) or die "cant open $annot_file...\n";
	open(OUT,">$chr.fa.all.GC") or die "cannot open $chr.fa.all.GC to write\n";
	print OUT "#Gene|Accession|Strand|Region\t#GC\tTotal\tRatio\n";
	while (<EX>) {
		chomp;
		my ($region,$rgType,$rgStart,$rgEnd) = split(/\t/);
		my $bases = substr($chr_seq, $rgStart, $rgEnd - $rgStart + 1);
		#print "$region,$rgType,$rgStart,$rgEnd: $bases\n";
		my $total = 1;
		if(defined($bases) &&  length($bases)>0) {
			$total = length($bases) ;
			my $gc = ($bases =~ s/[gc]//gi);
			my $ratio = sprintf("%.2f", $gc/$total);
			print OUT "$region\t$gc\t$total\t$ratio\n";
		}else{
			print OUT "$region\t0\t0\t0\n";
		}
	}
	close EX;
	close OUT;

}
