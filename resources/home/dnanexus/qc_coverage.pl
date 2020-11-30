#!/usr/bin/env perl 
#


use strict; 
use warnings; 

use TdtConfig; 

use Getopt::Long; 

use vars qw($opt_target $opt_file $opt_sample $opt_case); 
&GetOptions('target|t=s' => \$opt_target, 
	'file|f=s' => \$opt_file,
	'sample|s=s' => \$opt_sample,
	'case_control|c=s' => \$opt_case,
);

#print STDOUT "# Sample\tStatus\tCoverageValue\tCutoff\n"; 
my $conf = TdtConfig::readConfig('target', $opt_target); 
my $threshold = ($conf->{'COV_THRESHOLD'} ? $conf->{'COV_THRESHOLD'} : 20 );
my $level = ($threshold - 1);

open (my $fh, "<", $opt_file); 
my $header = <$fh>;
chomp $header;  
$header = "\t".$header;
my @header = split("\t", $header);  
my %header; 
for (my $i = 0; $i < scalar(@header); $i++){
	$header{$header[$i]} = $i; 
}

while (my $line = <$fh>){
	chomp $line; 
	next if ($line !~ m/^rcumsum/);
	my @line = split("\t", $line);
	my $val = $line[$header{$level}];
	my $cutoff =  $conf->{'MIN_COVGT'.$threshold.'X_'.$opt_case};
	my $status = ( $val > $cutoff) ? 'PASS' : 'FAIL';
	print STDERR "No cutoff defined\n" if (! defined $cutoff);  
	print STDOUT "$opt_sample\t$status\t$val\t$cutoff\n"; 
}
