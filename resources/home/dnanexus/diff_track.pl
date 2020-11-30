#!/usr/bin/env perl
#
# Calculates normalized coverage and coverage-difference.
#
# Params:
# 1 = Path to normal wiggle file (for the given chromosome)
# 2 = Path to tumor wiggle file (for the given chromosome)
# 3 = Average or median coverage of the normal sample
# 4 = Average or median coverage of the tumor sample
# 5 = Chromosome

use strict;
use warnings;

# Get params
my ($normal_wig, $tumor_wig, $normal_cov, $tumor_cov, $chr) = @ARGV; 

# Calculate coverage ratio
my $ratio = $normal_cov / $tumor_cov;

# Calculate normalized coverage and the coverage-difference and print
# these values to output files
open(normalized_fh,"> $chr-Dn.wig");
open(difference_fh,"> $chr-diff.wig"); 
open(tumor_wig_fh, $tumor_wig) or die "Error: Can't open $tumor_wig\n";
open(normal_wig_fh, $normal_wig) or die "Error: Can't open $normal_wig\n";
my ($normal_wig_cov, $tumor_wig_cov, $normalized_cov, $diff_cov);
while ($tumor_wig_cov = <tumor_wig_fh>) {
	chomp $tumor_wig_cov;
	$normal_wig_cov = <normal_wig_fh>;
	chomp $normal_wig_cov;
	# If this line doesn't contain coverage values, just print it out as it is
	if ($tumor_wig_cov =~ /^fixed/) {
		print normalized_fh "$tumor_wig_cov\n";
		print difference_fh "$tumor_wig_cov\n";
	}
	# Otherwise, if it's a coverage value, calculate normalized and diff coverage
	# and print them to output
	else {
		$normalized_cov = sprintf("%d", $tumor_wig_cov * $ratio);
		print normalized_fh "$normalized_cov\n";
		$diff_cov = $normalized_cov - $normal_wig_cov;
		print difference_fh "$diff_cov\n";
	}
}
close tumor_wig_fh;
close normal_wig_fh;
close normalized_fh;
close difference_fh;
