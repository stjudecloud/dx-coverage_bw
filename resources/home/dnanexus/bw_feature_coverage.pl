#!/bin/env perl
#
# a bit like bw_exon_coverage.pl, but broken down by feature
# within an isoform (intron/exon/utr)
#
# MNE 1/2015

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use File::Basename;
use Getopt::Long;
use List::Util qw(max);

use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;

use RefFlatFile();
use MiscUtils qw(dump_die median);
use FileUtils qw(read_simple_file);
# use DelimitedFile;
use Reporter;
use SJPreferredIsoform;
use TdtConfig;
use SampleName;
use DelimitedFile;

my %FLAGS;
my @GENES;

GetOptions(
	   \%FLAGS,

	   "-gene=s" => \@GENES,

	   "-genome=s",
	   "-refflat=s",
	   # override

	   "-preferred-isoforms=s",
	   # /nfs_exports/apps/gnu-apps/NextGen/nextgensupport2/NCBI/RefSeq/gene_transcript_matrix.withNM.mod
	   # add to genome config on research (there in clinical)

	   "-bams=s",
	   # BAMs to look up .bw files from

	   "-all-isoforms",
	   # don't restrict to SJ-preferred only


	   "-manifest=s",
	   # optional: check vs. Illumina manifest file

	  );

die "-gene" unless @GENES;

my $bw_files;
if (my $bl = $FLAGS{"bams"}) {
  $bw_files = bam_list_to_bw($bl);
}
die "need -bams [listfile] (or other method to get .bw list)" unless $bw_files;

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

my $f_pi = $FLAGS{"preferred-isoforms"} || die;
my $sjpi = new SJPreferredIsoform("-file" => $f_pi);

my %manifest;
if (my $mf = $FLAGS{manifest}) {
  # optional: map with Illumina ".bed"-like file
  my $df = new DelimitedFile("-file" => $mf,
			     "-headers" => 1,
			     );
  while (my $row = $df->get_hash()) {
    my $chr = $row->{Manifest_Chr} || die;
    push @{$manifest{$chr}}, $row;
  }
}

my %wanted_acc;
my %acc2gene;

foreach my $g (@GENES) {
  my $acc = $sjpi->get_preferred_isoform($g) || die;
  $wanted_acc{$acc} = 1;
  $acc2gene{$acc} = $g;
}

my $rf_fn = $FLAGS{refflat};
unless ($rf_fn) {
  my $genome = $FLAGS{genome} || die "specify -genome [genome] or -refflat [file]\n";
  die "specify -refflat (genome config lookup not implemented)";
}

my $rf = new RefFlatFile();
$rf->preserve_interbase(1);
$rf->parse_file(
		"-refflat" => $rf_fn,
		"-type" => "refgene",
	       );
my $rows_raw = $rf->rows();

my @wanted_rf;
my $isoform_setup_error;

if ($FLAGS{"all-isoforms"}) {
  #
  #  all isoforms for each gene
  #
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
} else {
  #
  # restrict to SJ-preferred isoform for each gene only
  #
  foreach my $r (@{$rows_raw}) {
    my $acc = $r->{name};
    if ($wanted_acc{$acc}) {
      push @wanted_rf, $r;
      $wanted_acc{$acc} = 2;
    }
  }

  foreach my $acc (sort keys %wanted_acc) {
    unless ($wanted_acc{$acc} == 2) {
      printf STDERR "ERROR: can't find %s (%s) in %s!\n", $acc, $acc2gene{$acc}, $rf_fn;
      $isoform_setup_error = 1;
    }
  }
}

die if $isoform_setup_error;


my @rpt_labels = qw(
		     bw
		     gene
		     accession
		     strand
		     chrom
		     start
		     end

		     feature
		     feature_number

		     coverage_median
		     coverage_max

		     coverage_percent_d20
		     coverage_percent_d30
		     coverage_percent_d40
		  );

if (%manifest) {
  push @rpt_labels, "manifest_overlap";
  push @rpt_labels, "manifest_contained";
}

foreach my $bw_file (@{$bw_files}) {
  printf STDERR "processing %s...\n", $bw_file;
  my $outfile = sprintf '%s.coverage.tab', basename($bw_file);
  my $bw = Bio::DB::BigFile->bigWigFileOpen($bw_file) || die;
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

      my $intervals = $bw->bigWigIntervalQuery($chrom, $start - 1 => $end);
      my @values;
      my $total_bases = 0;
      my $total_bases_less_20 = 0;
      my $total_bases_ge_20 = 0;
      my $total_bases_ge_30 = 0;
      my $total_bases_ge_40 = 0;

      for (my $i=$intervals->head;$i;$i=$i->next) {
	my $start = $i->start;
	my $end   = $i->end;
	my $val   = $i->value;
	#	    printf "%s\n", join " ", $start, $end, $val;
	$total_bases++;
	push @values, $val;
	if ($val < 20) {
#	  printf STDERR "only %d reads at %s %s.%d\n", $val, $fn, $chrom, $end;
	  $total_bases_less_20++;
	} else {
	  $total_bases_ge_20++;
	  $total_bases_ge_30++ if $val >= 30;
	  $total_bases_ge_40++ if $val >= 40;
	}
      }

      my %r;
      $r{bw} = basename($bw_file);
      $r{gene} = $rf_row->{gene} || die;
      $r{accession} = $rf_row->{name} || die;
      $r{strand} = $rf_row->{strand} || die;
      $r{feature} = $feature->{feature} || die;
      $r{feature_number} = $feature->{feature_number} || die;
      $r{chrom} = $chrom;
      $r{start} = $start;
      $r{end} = $end;

      $r{coverage_median} = median(\@values);
      $r{coverage_max} = @values ? max(@values) : 0;

#      my @nonzero = grep {$_ > 0} @values;
#      $r{coverage_nonzero_median} = @nonzero ? median(\@nonzero) : 0;


      $r{coverage_percent_d20} = sprintf '%4.1f', $total_bases_ge_20 * 100 / $total_bases;
      $r{coverage_percent_d30} = sprintf '%4.1f', $total_bases_ge_30 * 100 / $total_bases;
      $r{coverage_percent_d40} = sprintf '%4.1f', $total_bases_ge_40 * 100 / $total_bases;

      if (%manifest) {
	my $c = $chrom;
	$c =~ s/^chr//;
	my $rows = $manifest{$c} || die "no manifest entries for $c";
	my @manifest_overlap;
	my @manifest_contained;
	foreach my $mr (@{$rows}) {
	  my $ms = $mr->{Manifest_Start} || die;
	  my $me = $mr->{Manifest_End} || die;
	  my $m_label = sprintf '%s:%d-%d', $chrom, $ms, $me;

	  next if $start > $me;
	  next if $end < $ms;
	  # at this point we at least overlap
	  push @manifest_overlap, $m_label;
	  push @manifest_contained, $m_label if $start >= $ms and $end <= $me;
	}
	$r{manifest_overlap} = join ",", @manifest_overlap;
	$r{manifest_contained} = join ",", @manifest_contained;
      }

      $rpt->end_row(\%r);
    }
  }
  $rpt->finish();

}

sub bam_list_to_bw {
  my ($listf) = @_;

  my $list = read_simple_file($listf);
  my @bw;

  foreach my $bam (@{$list}) {
    my %info = SampleName::parse(basename($bam));

    my $bw = sprintf '/nfs_exports/genomes/1/projects/.allType/Eval-TSCA/.allTumor/SJ%s/NextGen/Amplicon/Coverage/%s/%s.bw',
      $info{disease},
	$info{subject},
	  $info{sample};
    push @bw, $bw;
  }
  return \@bw;
}
