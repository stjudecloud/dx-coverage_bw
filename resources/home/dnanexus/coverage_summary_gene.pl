#!/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use SampleName;
use Getopt::Long;
use TdtConfig;

#script to format results for coverage

my $result = GetOptions(
	"c|config=s" => \ my $config
);

my %conf = %{TdtConfig::readConfig(1,$config)};

my $sample=$conf{'SAMPLE'};
my $coverage_file=$conf{'COVERAGE_FILE'};
my $output_file=$conf{'OUTPUT_FILE'};
my $cutoff=$conf{'CUTOFF'};

if(!defined($cutoff)){
	$cutoff = 70;
}

my %parts = SampleName::parse($sample, "DIE");

my $disease = $parts{'disease'};
my $subject = $parts{'subject'};
my $type = $parts{'type'};


my %coverage;
open COVERAGE_FILE, $coverage_file or die("$coverage_file missing\n");
open OUTPUT, "> $output_file" or die("Cannot create output file $output_file");
print OUTPUT "Gene	Accession	Subject	Type	Coverage20x\n";
while(<COVERAGE_FILE>){
	chomp($_);
	my ($chr,$gene,$stand,$accession,$basesmorethan20x,$baseslessthan20x,$c20x,$c30x,$c40x,@exons) = split("\t",$_);
	if($chr ne "chr"){
		$c20x=~s/\%//g;
		$c30x=~s/\%//g;
		$c40x=~s/\%//g;
		$coverage{$gene}{$accession}{$subject}{$type}{'c20x'}=$c20x;
		$coverage{$gene}{$accession}{$subject}{$type}{'c30x'}=$c30x;
		$coverage{$gene}{$accession}{$subject}{$type}{'c40x'}=$c40x;
		if($c20x < $cutoff){
			print OUTPUT "$gene	$accession	$subject	$type	$c20x\n";
		}
	}
}
close COVERAGE_FILE;
close OUTPUT;

# this part was for producing a plot from R - probably not required
# print "done loading\n";
# print "list\tdonor\ttype\tgene\tcoverage_depth\tvalue\n";
# for my $gene (keys %coverage){
	# for my $donor (keys %{$coverage{$gene}}){
		# for my $type (keys %{$coverage{$gene}{$subject}}){
			# for my $coverage_type (keys %{$coverage{$gene}{$subject}{$type}}){
				# my $coverage=$coverage{$gene}{$subject}{$type}{$coverage_type};
				# $type=~s/[0-9]//g;
				# print $coverage_file."\t".$subject."\t".$type."\t".$gene."\t".$coverage_type."\t".$coverage."\n";
			# }
		# }
	# }
# }

