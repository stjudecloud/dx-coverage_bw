#!/usr/bin/env perl
###############################################################################
## Name:    exon_summary.pl                                                   #
## Author:  Gang Wu                                                           #
## Date:    March 21st, 2011                                                  #
## Purpose: get the mean coverage for all exons on all chromosomes            #
## This was based on exon_summary.pl, this version deals with *.exon.txt      #
## Usage:   perl exon_summary.pl F SJTALL004_G \
##   /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/chr2gene \
##  /nfs_exports/genomes/1/PCGP/BucketIntermediate/coverageWork/GRCh37/test/presum
###############################################################################

use strict;
use warnings;

#Check input parameters
die "need at least three arguments: gender, prefix, cov_dir\n" if(@ARGV<2);
my ($gender,$prefix,$annot_dir,$cov_dir) = @ARGV;
if ($ARGV[0] eq "M" || $ARGV[0] eq "F") {
	$gender = $ARGV[0];
}else{
	die "The gender can only be 'F' or 'M'!\n";
}


##get locations and GC content of each exons
print "get gene annotation..\n";
my %exons2gc;
my @chr2gc = glob("$annot_dir/*.all.GC");
foreach my $gc_file (@chr2gc) {
	open(G,"$gc_file") or die "cant open $gc_file\n";
	while (<G>) {
		chomp;
		if (!/^\#/) {
			my ($gae,$gc,$tot,$ratio) = split(/\t/);
			if ($gae =~ /exon\d+$/) {
				my ($gene,$accession,$strand,$seg) = split(/\|/,$gae);
				$seg =~ s/exon//;
				$exons2gc{$gae} = "$seg\t$ratio";
			}
		}
	}
}
my %gene2exons;
my @chr2genes = glob("$annot_dir/*.gene.txt");
foreach my $chr2gene (@chr2genes) {
	my ($chr) = ($chr2gene =~ /.+\/(chr.+)\.gene\.txt/);
	open(C,"$chr2gene") or die "cant open $chr2gene\n";
	while (<C>) {
		chomp;
		#A4GNT|NM_016161|-|UTR3  UTR     139325249       139325794
		#A4GNT|NM_016161|-|UTR5  UTR     139332789       139333919
		#A4GNT|NM_016161|-|exon1 exon    139333743       139333919) 
		my ($gae,$seg_type,$start,$end) = split(/\t/);
		if ($seg_type eq "Exon") {
			my ($gene,$accession,$strand,$seg) = split(/\|/,$gae);
			$gene2exons{$gae} = $gene."\t".$chr."\t".$strand."\t".$start."\t".$end;
			#print "$gae: $gene2exons{$gae}\n";
		}
	}
}

my %covs;
my $max_exon= 1;
my (%means,%chr2counts,%nr_means,%means_gc);
my @files = glob ("$cov_dir/$prefix*.Exon.txt");
foreach my $f (@files) {
	#grep ^# SJINF001_D_chr22_coverage.exon.txt
	##Gene|Accession|Strand|Segment 0 1 2 3 4 5 6 7 8 9 
	#A4GALT|NM_017436|-|exon1 48 27

	my ($pid,$chr) = ($f =~ /$cov_dir\/(.+)_(chr.+)_coverage/);
	
	print "processing $pid,$chr,$f...\n";
	open(F,$f) or die "cant open $f\n";
	while (<F>) {
		chomp;
		if (!/^\#/) {
			my ($gae,@counts) = split(/\s/);
			my ($gene,$accession,$strand,$exon) = split(/\|/,$gae);
			$exon =~ s/exon//;
			$max_exon = $exon if($exon>$max_exon);
			
			my ($total,$num);
			for (my $i=0; $i<=$#counts; $i++) {
				if ($counts[$i]>0) {
					if (defined($chr2counts{$chr}->{$i})) {
						$chr2counts{$chr}->{$i} += $counts[$i];
					}else{
						$chr2counts{$chr}->{$i} = $counts[$i];
					}
					$covs{$i} = 1;
				}
				$total += $counts[$i] * $i;
				$num   += $counts[$i];
			}

			my $mean = sprintf("%d",$total/$num);
			#print "$gene,$accession,$strand,$exon,$total,$num,$mean\n";
			$means{$chr."|".$gene."|".$strand."|".$accession}[$exon] = $mean; #
			$nr_means{$gene2exons{$gae}}= $mean; #
			$means_gc{$gae}=$mean;
		}
	}
}
 

print "output results..\n";

open(OUT,">$prefix.exon.coverage.means.txt") or die "cant open file to write\n";
foreach my $ga (sort keys %means) {
	print OUT "$ga";
	my @exon_means=@{$means{$ga}};
	for (my $exon=1; $exon<=$#exon_means; $exon++) {
		print OUT "\t$exon_means[$exon]";
	}
	print OUT "\n";
}
close OUT;


open(OUT,">$prefix.exon.coverage.counts.txt") or die "cant open file to write\n";
print OUT "Coverage";
my @chrs = (1..22,"X","Y");
@chrs = (1..22,"X") if($gender eq "F"); #Exclude Y chromosome if it is a female

foreach my $ch (@chrs) {
          print OUT "\tchr$ch";
}
print OUT "\tTotal\n";

foreach my $i (sort {$a<=>$b} keys %covs) {
	print OUT $i;
	my $total=0;
	foreach my $ch (@chrs) {
		if (defined($chr2counts{"chr".$ch}->{$i})) {
			print OUT "\t".$chr2counts{"chr".$ch}->{$i};
			$total += $chr2counts{"chr".$ch}->{$i};
		}else{
			print OUT "\t0";
		}
		
	}
	print OUT "\t$total\n";
}
close OUT;

open(OUT,">$prefix.exon.coverage.means.nr.sorted.txt") or die "cant open file to write\n";
foreach my $gae (sort {$nr_means{$a} <=> $nr_means{$b}}  keys %nr_means) {
	print OUT "$gae\t$nr_means{$gae}\n";
}
close OUT;

open(OUT,">$prefix.exon.coverage.means.with.GC.txt") or die "cant open file to write\n";
foreach my $gae (sort keys %means_gc) {
	if (defined($exons2gc{$gae})) {
		print OUT "$gae\t$exons2gc{$gae}\t$means_gc{$gae}\n";
	}
}
close OUT;
