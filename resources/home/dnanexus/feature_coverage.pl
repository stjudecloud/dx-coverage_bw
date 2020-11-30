#!/bin/env perl
#
# report isoform features and coverage provided by an external source
# (e.g. bed).  See also bw_feature_coverage.pl.
#
# MNE 3/2015

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;

use MiscUtils qw(dump_die);
use RefFlatFile;
# use DelimitedFile;
# use Reporter;

my %FLAGS;
#GetOptions(\%FLAGS, );

my @GENES;
my @BED;
GetOptions(
	   \%FLAGS,
	   "-gene=s" => \@GENES,

	   "-genome=s",
	   "-refflat=s",
	   # override

	   "-all-isoforms",
	   # don't restrict to SJ-preferred only

	   # external coverage source:
	   "-bed=s" => \@BED,
	   "-bed-glob",
	  );

die "-gene" unless @GENES;
push @BED, glob("*.bed") if $FLAGS{"bed-glob"};
die "-bed [file] or -bed-glob" unless @BED;

my $rf_fn = $FLAGS{refflat} || die "-refflat";
# e.g.
# /nfs_exports/genomes/1/Homo_sapiens/Hg19/mRNA/RefSeq/refFlat.15AUG2014.txt

my $rf = new RefFlatFile();
$rf->preserve_interbase(1);
$rf->parse_file(
		"-refflat" => $rf_fn,
		"-type" => "refgene",
	       );
my $rows_raw = $rf->rows();

my @wanted_rf;
my $isoform_setup_error;

my %wanted_genes = map {$_, 1} @GENES;
foreach my $r (@{$rows_raw}) {
  my $gene = $r->{gene} || die;
  if ($wanted_genes{$gene}) {
    push @wanted_rf, $r;
    $wanted_genes{$gene} = 2;
  }
}
foreach my $g (sort keys %wanted_genes) {
  unless ($wanted_genes{$g} == 2) {
    printf STDERR "ERROR: no isoforms found for %s\n", $g;
    $isoform_setup_error = 1;
  }
}

my $coverage = parse_covered_from_bed();

my @rpt_labels = qw(
		     gene
		     accession
		     strand
		     chrom
		     start
		     end
		     feature
		     feature_number

		     coverage_nt
		     coverage_percent
		     coverage_touched
		  );

my $outfile = "coverage.tab";

my $rpt = new Reporter(
		       "-file" => $outfile,
		       "-delimiter" => "\t",
		       "-labels" => \@rpt_labels,
		       "-auto_qc" => 1,
		      );

foreach my $rf_row (@wanted_rf) {
  my $features = $rf->get_feature_summary("-row" => $rf_row);
  my $chrom = $rf_row->{chrom} || die;

  foreach my $feature (@{$features}) {
    my $start = $feature->{start};
    my $end = $feature->{end};
    # 1-based

    my %r;
    $r{gene} = $rf_row->{gene} || die;
    $r{accession} = $rf_row->{name} || die;
    $r{strand} = $rf_row->{strand} || die;
    $r{feature} = $feature->{feature} || die;
    $r{feature_number} = $feature->{feature_number} || die;
    $r{chrom} = $chrom;
    $r{start} = $start;
    $r{end} = $end;

    my %covered;
    my $cr = $coverage->{$chrom} || die;
    my %touched;
    for (my $i = $start; $i <= $end; $i++) {
      my $set = $cr->{$i};
      if ($set) {
	$covered{$i} = 1;
	foreach my $entry (@{$set}) {
	  my $touched = sprintf '%s:%d-%d', @{$entry};
	  $touched{$touched} = 1;
	}
      }
    }

    my $ccount = scalar keys %covered;
    my $len = $end - $start + 1;
    my $cpct = sprintf "%.2f", $ccount * 100 / $len;
    $r{coverage_nt} = $ccount;
    $r{coverage_percent} = $cpct;
    $r{coverage_touched} = join ",", sort keys %touched;

    $rpt->end_row(\%r);
  }
}
$rpt->finish();




sub parse_covered_from_bed {
  die unless @BED;
  my %coverage;
  foreach my $bed (@BED) {
    printf STDERR "parse %s\n", $bed;
    open(IN, $bed) || die "can't open $bed";
    while (<IN>) {
      chomp;
      my @f = split /\t/, $_;
      next unless @f > 1;
      die unless @f == 4;
      my ($chr, $start, $end, $thing) = @f;
      $start++;
      # convert from interbase to in-base
      die unless $chr =~ /^chr/;

      my $ref = [ $chr, $start, $end ];

      for (my $i = $start; $i <= $end; $i++) {
	push @{$coverage{$chr}{$i}}, $ref;
      }
    }
    close IN;
  }
  return \%coverage;
}

