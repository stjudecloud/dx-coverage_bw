#!/usr/bin/perl
###############################################################################
## Name:    chr2gene.pl                                                       #
## Author:  Gang Wu                                                           #
## Date:    Dec. 22, 2010                                                     #
## Purpose: Split the refFlat file to gene entity annotation file per chr     #
##          This script includes all UTRs, introns and exons                  #
###############################################################################

use strict;
use warnings;

if  (@ARGV<1) {
   print STDERR "Usage: chr2gene.pl <refFlat.txt>\n"; 
   exit;
}
my $refFlat = $ARGV[0];
#"../mm9.refFlat.Dec.22.2010.txt"; #"/nfs_exports/apps/gnu-apps/NextGen/nextgensupport/refFlat.ucsc.txt";
##geneName       name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds
#FAM138F NR_026820       chr1    -       24474   25944   25944   25944   3       24474,25139,25583,      25037,25344,25944,
my %chr2gene;
my %accessions;

open(F,$refFlat) or die "cant open $refFlat\n";
while (<F>) {
	chomp;
	if (!/^\#/) { #skip header	
		my ($gene,$accession,$chr,$strand,$txS,$txE,$cdsS,$cdsE,$exonCount,$exonStarts,$exonEnds) = split(/\t/);
		my $acc = $accession;
		$acc = "$acc.$accessions{$accession}" if (defined($accessions{$accession}));
		$accessions{$accession}++;
		my $gas = $gene."|".$acc."|".$strand;
		next if($chr =~ /_/);
		if ($strand eq "+") {
			if ($txS<$cdsS) { #5'-UTR
				$chr2gene{$chr}->{$gas."|UTR5\tUTR"} = "$txS\t".($cdsS-1);
			}
			if ($cdsE<$txE) { #3'-UTR
				$chr2gene{$chr}->{$gas."|UTR3\tUTR"} = ($cdsE+1)."\t$txE";
			}
		} elsif($strand eq "-"){
			if ($txS<$cdsS) { #3'-UTR
				$chr2gene{$chr}->{$gas."|UTR3\tUTR"} = "$txS\t".($cdsS-1);
			}
			if ($cdsE<$txE) { #5'-UTR
				$chr2gene{$chr}->{$gas."|UTR5\tUTR"} = ($cdsE+1)."\t$txE";
			}
		} else{
			die "no strand info for $accession!\n";
		}

		##++++++++exons and introns++++++++
		chop($exonStarts);
		my @exStarts = split(/,/,$exonStarts);
		chop($exonEnds);
		my @exEnds   = split(/,/,$exonEnds);
		# the following method won't work if the CDS starts in the 2nd or even 3rd exons. The same problem for CDS ends
		#$exStarts[0] = $cdsS; #replace the start of the first exon with CDS start
		#$exEnds[$#exEnds] = $cdsE; #replace the end of the last exon with CDS end
		for (my $i=0; $i<=$#exStarts; $i++) {
			if($strand eq "+"){ #forward strand
				my $exon = $gas."|exon".($i+1)."\tExon";
				$chr2gene{$chr}->{$exon} = "$exStarts[$i]\t$exEnds[$i]";
				if ($i>0) { # at least two exons
				  my $intron = $gas."|intron".$i."\tIntron";
				  $chr2gene{$chr}->{$intron} = ($exEnds[$i-1]+1)."\t".($exStarts[$i]-1) if($exEnds[$i-1] < $exStarts[$i]);
				}

				#for coding exons
				# the assumption is cdsS<cdsE AND $exStarts[$i] < $exEnds[$i]
				#case 1: left out (meaning the current exon is completely outside of CDS and on the left
				if ($exEnds[$i] <= $cdsS) {
					#skip; this includes a single base within CDS, which is ignored
				} elsif ($exStarts[$i] < $cdsS && $exEnds[$i] < $cdsE){ 
				#case 2: left overlaps, partially overlapping at the left side of CDS
					my $cdsexon = $gas."|exon".($i+1).".5partial\tCodingExon";
					my $s = $cdsS;
					my $e = $exEnds[$i];
					$chr2gene{$chr}->{$cdsexon} = "$s\t$e";
				}
				elsif ($exStarts[$i] >= $cdsS && $exEnds[$i] <= $cdsE){
				#case 3: inside, completely in the CDS, for majority of exons, extreme case=>completely overlapping
					my $cdsexon = $gas."|exon".($i+1).".within\tCodingExon";
					$chr2gene{$chr}->{$cdsexon} = "$exStarts[$i]\t$exEnds[$i]";
				}
				elsif ($exStarts[$i] > $cdsS && $exStarts[$i] < $cdsE && $exEnds[$i] > $cdsE){
				#case 4: right overlaps, partially overlapping at the right side of CDS
					my $cdsexon = $gas."|exon".($i+1).".3partial\tCodingExon";
					my $s = $exStarts[$i];
					my $e = $cdsE;
					$chr2gene{$chr}->{$cdsexon} = "$s\t$e";
				}
				elsif ($exStarts[$i] >= $cdsE){
				#case 5: right out (meaning the current exon is completely outside of CDS and on the right side
					#skip
				}
				elsif ($exStarts[$i] < $cdsS && $exEnds[$i] > $cdsE){
				#case 6: inclusive, the whole CDS is within the current exon, for the case of long exon and short CDS
					my $cdsexon = $gas."|exon".($i+1).".inclusive\tCodingExon";
					$chr2gene{$chr}->{$cdsexon} = "$cdsS\t$cdsE";
				}
				else{
					#skip
				}

			}elsif($strand eq "-"){ #reverse strand
				my $exon = $gas."|exon".($#exStarts - $i + 1)."\tExon";
				$chr2gene{$chr}->{$exon} = "$exStarts[$i]\t$exEnds[$i]";
				if ($i>0) { # at least two exons
				  my $intron = $gas."|intron".($#exStarts - $i + 1)."\tIntron";
				  $chr2gene{$chr}->{$intron} = ($exEnds[$i-1]+1)."\t".($exStarts[$i]-1) if($exEnds[$i-1] < $exStarts[$i]);
				}

				#for coding exons
				# the assumption is cdsS<cdsE AND $exStarts[$i] < $exEnds[$i]
				#case 1: left out (meaning the current exon is completely outside of CDS and on the left
				if ($exEnds[$i] <= $cdsS) {
					#skip; this includes a single base within CDS, which is ignored
				} elsif ($exStarts[$i] < $cdsS && $exEnds[$i] < $cdsE){ 
				#case 2: left overlaps, partially overlapping at the left side of CDS
					my $cdsexon = $gas."|exon".($#exStarts - $i + 1).".3partial\tCodingExon";
					my $s = $cdsS;
					my $e = $exEnds[$i];
					$chr2gene{$chr}->{$cdsexon} = "$s\t$e";
				}
				elsif ($exStarts[$i] >= $cdsS && $exEnds[$i] <= $cdsE){
				#case 3: inside, completely in the CDS, for majority of exons, extreme case=>completely overlapping
					my $cdsexon = $gas."|exon".($#exStarts - $i + 1).".within\tCodingExon";
					$chr2gene{$chr}->{$cdsexon} = "$exStarts[$i]\t$exEnds[$i]";
				}
				elsif ($exStarts[$i] > $cdsS && $exStarts[$i] < $cdsE && $exEnds[$i] > $cdsE){
				#case 4: right overlaps, partially overlapping at the right side of CDS
					my $cdsexon = $gas."|exon".($#exStarts - $i + 1).".5partial\tCodingExon";
					my $s = $exStarts[$i];
					my $e = $cdsE;
					$chr2gene{$chr}->{$cdsexon} = "$s\t$e";
				}
				elsif ($exStarts[$i] >= $cdsE){
				#case 5: right out (meaning the current exon is completely outside of CDS and on the right side
					#skip
				}
				elsif ($exStarts[$i] < $cdsS && $exEnds[$i] > $cdsE){
				#case 6: inclusive, the whole CDS is within the current exon, for the case of long exon and short CDS
					my $cdsexon = $gas."|exon".($#exStarts - $i + 1).".inclusive\tCodingExon";
					$chr2gene{$chr}->{$cdsexon} = "$cdsS\t$cdsE";
				}
				else{
					#skip
				}

			}else{
				die "no strand info for $accession!..\n";
			}
		}
	}
}
close F;


#output
foreach my $chr (sort keys %chr2gene) {
	open(OUT,">$chr.gene.txt") or die "cant open file to write\n";
	my %segments = %{$chr2gene{$chr}};
	foreach my $seg (sort keys %segments) { #sort numerically
		print OUT "$seg\t$segments{$seg}\n";
	}
	close OUT;
}

