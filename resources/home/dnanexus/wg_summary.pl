#!/usr/bin/env perl
###############################################################################
## Name:    wg_summary.pl                                                     #
## Author:  Gang Wu                                                           #
## Date:    July 29th, 2010                                                   #
## Purpose: Combine the coverage summary across chromosomes                   #
## Changes: This version no longer needs paired germline and diagnosis samples#
## Usage:   perl wg_summary.pl M SJMB002_D_102927016 .                        #
###############################################################################

use strict;
use warnings;
my (%covs,%chr2counts);

#Check input parameters
die "need at least three arguments: gender, prefix, cov_dir\n" if(@ARGV<2);
my ($gender,$prefix,$cov_dir) = @ARGV;
if ($ARGV[0] eq "M" || $ARGV[0] eq "F") {
	$gender = $ARGV[0];
}else{
	die "The gender can only be 'F' or 'M'!\n";
}

#Process coverage summary files
my @files = glob ("$cov_dir/$prefix*.WG.txt");
foreach my $f (@files) {
	#head SJINF001_D_chr4_coverage_counts_wg.txt
	#Coverage        Sites
	##N      3976000
	#0       218358
	my ($pid,$chr) = ($f =~ /$cov_dir\/(.+)_(chr.+)_coverage/);
	
	print "processing $pid,$chr,$f...\n";
	open(F,$f) or die "cant open $f\n";
	while (<F>) {
		chomp;
		if (/^\d/) {
			my ($cov,$count) = split(/\t/);
			$chr2counts{$chr}->{$cov} = $count;
			$covs{$cov} = 1;
		}
	}
}
 

print "output results..\n";

open(OUT,">$prefix.wg.coverage.counts.txt") or die "cant open file to write\n";
print OUT "Coverage";
my @chrs = (1..22,"X","Y");
@chrs = (1..22,"X") if($gender eq "F"); #Exclude Y chromosome if it is a female
foreach my $ch (@chrs) {
	
	print OUT "\tchr$ch";

}
print OUT "\tTotal\n";

foreach my $cov (sort {$a<=>$b} keys %covs) {
	print OUT $cov;
	my $total=0;
	foreach my $ch (@chrs) {
			if (defined($chr2counts{"chr".$ch}->{$cov})) {
				print OUT "\t".$chr2counts{"chr".$ch}->{$cov};
				$total += $chr2counts{"chr".$ch}->{$cov};
			}else{
				print OUT "\t0";
			}
	}	
	print OUT "\t$total\n";
}
close OUT;
