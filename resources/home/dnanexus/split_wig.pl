use strict;
use warnings;
my$current_chr;

# Input file from command line
my $infile = $ARGV[0];
open ( INFILE, "< $infile" ) or die "Can't open $infile : $!";

# Initial output file
my$outfile_stem=0;
my$outfile="$infile."."$outfile_stem".".wig";
open( OUTFILE, "> $outfile" ) or die "Can't open $outfile : $!";

while(<INFILE>){

	chomp;
	if ( $_ =~ m/step/){
	
		$current_chr=$_;
		close OUTFILE ;		
		print "$current_chr\n";
		$outfile_stem++;
		$outfile="$infile."."$outfile_stem".".wig";
		print "$outfile\n";
		open( OUTFILE, "> $outfile" ) or die "Can't open $outfile : $!";
		print OUTFILE "$_\n"; 
		
		}
	
	else{			
		print OUTFILE "$_\n";
		}
	
	}

close INFILE;	
close OUTFILE;	