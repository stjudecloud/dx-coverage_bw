#!/usr/bin/env perl
###############################################################################
## Name:    bw.covsum.pl                                                      #
## Author:  Gang Wu                                                           #
## Date:    Dec. 22, 2010                                                     #
## Purpose: Count the number of sites of each coverage for NGS data           #
##           using UCSC refFlat for gene annotations                          #
## Usage:    perl bw_covsum.pl -chr 22 \
##      -cov_file /nfs_exports/genomes/1/PCGP_dev/HematopoieticMalignancies/SJINF/NextGen/WholeGenome/Coverage/SJINF001/SJINF001_D.bw \
##      -annot_dir /nfs_exports/apps/gnu-apps/NextGen/nextgensupport/WashU_hg18_nib \
##      -chr_size /nfs_exports/apps/gnu-apps/NextGen/nextgensupport/hg18/hg18.sizes \
##      -out_dir . \                                                          #
##      -sample SJINF001_D \                                                  #
##      -seg_type WG                                                          #
## Assuming the coverge file is named as 'sample.bw'                          #
## Revision: to use UCSC bigWiggle file.                                      #
###############################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;

my %opt;
GetOptions(\%opt, 'cov_file=s',
                  'annot_dir=s',
                  'chr_size=s',
                  'chr=s',
                  'sample=s',
                  'out_dir=s',
                  'seg_type=s',
                  'help|?');
usage() if $opt{'help'} || keys(%opt) < 1 ;

#################################################################
####                        MAIN
##################################################################

####+++++set the global parameters ++++++++
my ($chr,$chr_index,$coverage_file,$annotation_file_path,$sample,$out_dir,$seg_type,%chr_sizes);
set_parms();
# hg18 chromosome lengths #plan to save as a parameter file
#my %chr_sizes = #(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, #1-6
#		 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,         #7-12
#		 114142980, 106368585, 100338915, 88827254, 78774742,76117153,             #13-18
#		 63811651, 62435964, 46944323, 49691432, 154913754, 57772954, 16571);      #19-25


####+++++set the annotation of genes for each base++++++++
#my @seg_types = ("Exon","UTR","CodingExon","Intron","WG");
#"Intron" takes too much memory, suggest to use a dedicated node to process it. (-pe smp 7)
#each segment can be a utr, an segment, an intron or a coding segment
if($seg_type eq "WG"){
	my $b2s_ref = set_wg_annotation($seg_type);
	process_genome($seg_type,$b2s_ref,$sample);
}else{
#	my $b2s_ref = set_gene_annotation($seg_type);
	process_gene($seg_type,$sample);
}


#################### Subroutines ######################
#######################################################
# Title   : usage()
# Function: message about the program and how to use it
#######################################################
sub usage{
    print << "EOF";

    This program summarizes NGS coverage by taking the chromosome number, 
    patient ID, diretory of coverage files and annotation files

    usage: $0 [-help] [options]

     -help         : this (help) message
     -chr          : the chromosome to process, a sing value in: 1..22, X, Y
     -chr_size     : the file with two columns: "chr1", chr_length
     -cov_file     : the UCSC bigWiggle file
     -annot_dir    : the directory of annotation files
     -sample       : the name of the patient
     -out_dir      : the directory to write results
     -seg_type     : "Exon","CodingExon","UTR", "Intron", "WG", or "Promoter"

EOF
    exit;
}
#######################################################
# Title   : set_parms
# Function: set global parameters
#######################################################
sub set_parms{
	my %chr_list;
	#check chr_size value
	if (!defined($opt{'chr_size'})){
		die "-chr_size undefined!\n";
	}elsif(!(-e $opt{'chr_size'})){
		die "chr_size file does not exist!\n";
	}else{
		open(L,$opt{'chr_size'}) or die "cant open $opt{'chr_size'}\n";
		while (<L>) {
			chomp;
			my ($c,$s) = split;
			$chr_sizes{$c} = $s;
			$c =~ s/chr//;
			$chr_list{$c} =  "valid_chromosome";
		}
	}
	
	#check chr value
	if (!defined($opt{'chr'})){
		die "-chr undefined!\n";
	}elsif(!defined($chr_list{$opt{'chr'}})){
		die "wrong value for -chr! \nAccepted values: 1..22, X, Y, M\n";
	}else{
		$chr = $opt{'chr'};
	}

	my %seg_types = ("Exon"  => 1, "CodingExon" => 1,
			"UTR"    => 1, "Intron"     => 1,
			"WG"     => 1, "Promoter"   => 1);

	#check seg_type  value
	if (!defined($opt{'seg_type'})){
		die "-seg_type undefined!\n";
	}elsif(!defined($seg_types{$opt{'seg_type'}})){
		die "wrong value for -seg_type! \nAccepted values: Exon, Intron, CodingExon, UTR, WG, Promoter\n";
	}else{
		$seg_type = $opt{'seg_type'};
	}

	#check paths
	if (!defined($opt{'cov_file'})){
		die "-cov_file undefined!\n";
	}elsif (-e $opt{'cov_file'}) {
		$coverage_file = $opt{'cov_file'};
	}else{
		die "cov_file $opt{'cov_file'} does not exist!\n";
	}

	if (!defined($opt{'annot_dir'})){
		die "-annot_dir undefined!\n";
	}elsif (-e $opt{'annot_dir'}) {
		$annotation_file_path = $opt{'annot_dir'};
	}else{
		die "annot_dir $opt{'annot_dir'} does not exist!\n";
	}

	if (!defined($opt{'sample'})){
		die "-sample undefined!\n";
	}elsif (-e $opt{'cov_file'}) {
		$sample = $opt{'sample'};
	}else{
		die "coverage file directory ".$opt{'cov_file'}." does not exist!\n";
	}

	if (!defined($opt{'out_dir'})){
		die "-out_dir undefined!\n";
	}elsif (-e $opt{'out_dir'}) {
		$out_dir = $opt{'out_dir'};
	}else{
		die "out_dir $opt{'out_dir'} does not exist!\n";
	}
	
	$annotation_file_path =~ s/\/$//; #remove the ending forward slash if exists.
	$out_dir =~ s/\/$//; #remove the ending forward slash if exists.
}

#######################################################
# Title   : set_gene_annotation
# Function: set the annotation for each genic base
#######################################################
sub set_gene_annotation{
	my $seg_type = shift;
	my @base2segments;
	my $bases = 0;
	my $gene_file = "$annotation_file_path/chr$chr.gene.txt";
	print "  set $seg_type annotation for each base on chr$chr with $gene_file\n";
	open(EX, $gene_file) or die "cant open $gene_file...\n";
	while (<EX>) {
		chomp;
		my ($segment,$stype,$segStart,$segEnd) = split(/\t/);
		if($stype eq $seg_type){ #only process a particular segment type per run
			for (my $base=$segStart; $base<=$segEnd; $base++) { #start and end are inclusive
				if (defined($base2segments[$base])) { #when a base belongs to more than 2 segments in the same or different genes
					$base2segments[$base] = $base2segments[$base].",".$segment;
					$bases++;
				}else{
					$base2segments[$base] = $segment;
					$bases++;
				}
			}
		}
	}
	close EX;
	
	print "  number of bases in $seg_type on chr$chr: $bases\n";
	return \@base2segments;
}

#######################################################
# Title   : set_wg_annotation
# Function: get the position of Ns (gaps) in the whole genome
#######################################################
sub set_wg_annotation{
	my $seg_type = shift;
	my %base2ns;
	my $bases = 0;
	my $ns_file = "$annotation_file_path/$chr.fa.Ns";
	print "  get sequencing gaps on chr$chr with $ns_file\n";
	open(NS, $ns_file) or die "cant open $ns_file...\n";
	while (<NS>) {
		chomp;
		my ($start,$end)=split(/,/);
		for (my $i=$start; $i<=$end; $i++) {
			$base2ns{$i} = 1;
		}
	}
	close NS;
	my @ns = keys %base2ns;
	print "  number of Ns on chr$chr: ". ($#ns + 1) . "\n";
	return \%base2ns;
}


#######################################################
# Title   : process_genome
# Function: process the coverage files for whole genome
#######################################################
sub process_genome {
	my $seg_type=shift;
	my $b2s_ref=shift;
	my $prefix = shift;
	my %base2ns = %{$b2s_ref}; #de-referencing the array reference
	my @ns = keys %base2ns;
	
	print "  read coverage from $coverage_file\n";
	my ($total,$informative, %counts);
	my $max_cov = 1;
	my $wig = Bio::DB::BigFile->bigWigFileOpen($coverage_file) || die;
	# query each of the intervals (fixed or variable step values)
	my $intervals = $wig->bigWigIntervalQuery("chr$chr", 0 => $chr_sizes{"chr$chr"});
	for (my $i=$intervals->head;$i;$i=$i->next) {
		my $val         = $i->value;
		my $base_number = $i->end;
		if (!defined($base2ns{$base_number})) { #excluding sequencing gaps
			my $cov = $val;
			$informative++;
			$counts{$cov}++;
			$max_cov = $cov if($cov > $max_cov);
		}
		$total++;
	}
	print "\ttotal: $total\tinformative: $informative\n";

	my $out_file = "$out_dir/$prefix"."_chr$chr"."_coverage.$seg_type.txt";
	open(OUT,">$out_file") or die "cant open $out_file to write!\n";
	print OUT "Coverage\tSites\n";		
	print OUT "#N\t".($#ns + 1) . "\n";
	foreach my $cov (sort {$a <=> $b} keys %counts) {
		print OUT "$cov\t$counts{$cov}\n";
	}
	close OUT;

}


#######################################################
# Title   : process_gene
# Function: process the coverage files for each segment type
#######################################################
sub process_gene {
	my $seg_type=shift;
	my $prefix = shift;

	print "  output results $prefix\n";
	my $out_file = "$out_dir/$prefix"."_chr$chr"."_coverage.$seg_type.txt";
	open(OUT,">$out_file") or die "cant open $out_file to write!\n";
	print OUT "Segment";
	for (my $i=0; $i<=100; $i++){
		print OUT " $i";
	}
	print OUT "\n";
	print "  read coverage from $coverage_file\n";
	my $wig = Bio::DB::BigFile->bigWigFileOpen($coverage_file) || die;
	# query each of the intervals (fixed or variable step values)
	my $gene_file = "$annotation_file_path/chr$chr.gene.txt";
	print "  get $seg_type annotation for each base on chr$chr with $gene_file\n";
	open(EX, $gene_file) or die "cant open $gene_file...\n";
	while (<EX>) {
		chomp;
		my ($segment_max,@counts);
		my ($segment,$stype,$segStart,$segEnd) = split(/\t/);
		if($stype eq $seg_type){ #only process a particular segment type per run
			my $intervals = $wig->bigWigIntervalQuery("chr$chr", ($segStart-1) => $segEnd);
			for (my $i=$intervals->head;$i;$i=$i->next) {
				my $cov = $i->value;
				if (defined($segment_max)) {
					$segment_max = $cov if($segment_max < $cov);
				}
				else{
					$segment_max = $cov;
				}
				$counts[$cov]++;
			}
				
			print OUT $segment;
			for (my $i=0; $i<=$segment_max; $i++) {
				if (defined($counts[$i])) {
					print OUT " $counts[$i]";
				}else{
					print OUT " 0";
				}
			}
			print OUT "\n";
		}
	}
	close EX;
	close OUT;
}
