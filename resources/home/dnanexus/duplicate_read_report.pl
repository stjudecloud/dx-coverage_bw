#!/bin/env perl
# find reads with high repeated read coverage
#
# limited read CTCAGAGAGAACCCACCATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCG at 16 pos 222895 with count 739961
# /rgs01/resgen/prod/tartan/index/data/PanTARGET/PanTARGET/SJAML040602_D1/TRANSCRIPTOME/bam/SJAML040602_D1.bam
#
# MNE 11/2015

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use Bio::DB::Sam;
use Bio::DB::Sam::Constants;

use Cluster;
use CommandLineRebuilder;
use ReferenceNameMapper;
use MiscUtils qw(dump_die build_argv_list);
use Reporter;
use GeneAnnotation;
use DelimitedFile;
use SampleTumorNormal;

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-bam=s",
	      "-limit=i",

	      "-out=s",

	      "-ref=s",
	      "-start=s",
	      "-end=s",
	      # optional restrict

	      "-ping=i",
	      "-sublist=i",

	      "-force",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if (my $size = $FLAGS{sublist}) {
  generate_sublist($size);
  exit(0);
}
my $LIMIT = $FLAGS{limit} || die "-limit";
my $genome = $FLAGS{genome} || die "-genome";
my $ping = $FLAGS{ping};

my @f_bam;
if (my $b = $FLAGS{bam}) {
  push @f_bam, $b;
}
die unless @f_bam;

my $GA = new GeneAnnotation(
#			    "-style" => "refgene",
			    "-style" => "refgene_flatfile",
			    "-genome" => $genome,
			   );

my $READ_COUNT;
my $LAST_CHR;
my $LAST_START;
my $OUT_COUNT;
my %TRACK;
my $RPT;

foreach my $f_bam (@f_bam) {
  my ($bam, $rnm) = bam_setup($f_bam);

  my ($chr, $start, $end);
  if ($start = $FLAGS{start}) {
    $end = $FLAGS{end} || $start;
  }
  if (my $chr_raw = $FLAGS{ref}) {
    $chr = $rnm->find_name($chr_raw) || die "no BAM chrom for $chr_raw";
  }

  my $interval = "";
  if ($chr) {
    $interval = $chr;
    if ($start) {
      $interval .= sprintf ":%d-%d", $start, $end;
    }
  }

  my $outfile = $FLAGS{out};
  unless ($outfile) {
    $outfile = sprintf "%s", basename($f_bam);
    $outfile .= "." . $chr if $chr;
    $outfile .= ".duplicate_report.tab", basename($f_bam);
  }

  next if -s $outfile and not($FLAGS{force});

  $RPT = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   chr
					   start
					   sequence
					   count
					   gene
					)
				      ],
			 "-auto_qc" => 1,
			);

  $READ_COUNT = 0;
  $OUT_COUNT = 0;
  $LAST_CHR = "";
  $LAST_START = -1;
  %TRACK = ();

  if ($interval) {
    $bam->fetch($interval, \&process_one);
    # the -interval parameter does not seem to work properly for
    # large intervals (e.g. whole chromosomes) in the installed version
    # of bioperl.  "fetch" seems to work OK, though requires a callback sub
    # to process each read  :/
  } else {
    my $iterator = $bam->features("-iterator" => 1);
    while (my $align = $iterator->next_seq()) {
      process_one($align);
    }
  }

  flush();
  $RPT->finish();
}


sub process_one {
  # ugh: since this is a callback called by third-party code,
  # we have to use global variables
  my ($a) = @_;
  $READ_COUNT++;
  if ($ping and $READ_COUNT % $ping == 0) {
    printf STDERR "%s: processed %d at %s.%d, detected=%d\n", scalar(localtime), $READ_COUNT, $LAST_CHR, $LAST_START, $OUT_COUNT;
  }
  if (my $this_start = $a->start()) {
    my $this_chr = $a->seq_id();
    if ($a->get_tag_values("DUPLICATE")) {
      flush() if $this_chr ne $LAST_CHR or $this_start != $LAST_START;
      my $seq = $a->query->dna();
      $TRACK{$seq}++;
    }
    $LAST_CHR = $this_chr;
    $LAST_START = $this_start;
  }
}

sub flush {
  foreach my $seq (keys %TRACK) {
    my $count = $TRACK{$seq};
    if ($count > $LIMIT) {
      $OUT_COUNT++;
      my $gene = "";
      if ($GA->find(
		    "-reference" => $LAST_CHR,
		    "-start" => $LAST_START,
		    "-end" => $LAST_START
		   )) {
	my %genes;
	foreach my $h (@{$GA->results_rows}) {
	  my $g = $h->{name2} || die;
	  $genes{$g} = 1 if $g;
	}
	$gene = join ",", sort keys %genes;
      }

      my %r;
      $r{chr} = $LAST_CHR;
      $r{start} = $LAST_START + 1;
      # convert to 1-based
      $r{sequence} = $seq;
      $r{count} = $count;
      $r{gene} = $gene;

#      dump_die(\%r, "Debug", 1);

      $RPT->end_row(\%r);
    }
  }

  %TRACK = ();
}

sub bam_setup {
  my ($f_bam) = @_;
  die "where is $f_bam" unless -s $f_bam;

  my $bam = Bio::DB::Sam->new(
			      "-bam" => $f_bam,
			      "-expand_flags" => 1,
			     );
  my $rnm = new ReferenceNameMapper();
  foreach my $chrom ($bam->seq_ids()) {
    # standardize reference names
    $rnm->add_name($chrom);
  }
  return ($bam, $rnm);
}

sub generate_sublist {
  my ($size) = @_;
  my @reports = glob("*duplicate_report.tab");
  die unless @reports;
  my $rpt;
  my %all;
  my $df;
  foreach my $r (@reports) {
    if (-s $r) {
      $df = new DelimitedFile("-file" => $r,
			      "-headers" => 1,
			     );
      unless ($rpt) {
	my $outfile = sprintf "sublist_%d.tab", $size;
	$rpt = $df->get_reporter(
				 "-file" => $outfile,
				 "-extra" => [ "disease", "sample" ],
				 "-extra-prepend" => 1,
				);
      }
      while (my $row = $df->get_hash()) {
	if ($row->{count} >= $size) {
	  my $stn = new SampleTumorNormal();
	  $stn->parse($r) || die;

	  my $sample = basename($r);
	  $sample =~ s/\.bam.*//;
	  $row->{sample} = $sample;
	  $row->{disease} = $stn->disease;

	  $rpt->end_row($row);

	  my $key = $row->{gene} || $row->{chr} || die;
	  push @{$all{$sample}{$key}}, $row;
	}
      }
    }
  }
  $rpt->finish();

  my $outfile = sprintf "sublist_%d_single.tab", $size;
  $rpt = $df->get_reporter(
			   "-file" => $outfile,
			   "-extra" => [ "disease", "sample" ],
			   "-extra-prepend" => 1,
			  );
  foreach my $sample (sort keys %all) {
    foreach my $key (sort keys %{$all{$sample}}) {
      my $rows = $all{$sample}{$key};
      my @sorted = sort {$b->{count} <=> $a->{count}} @{$rows};
      $rpt->end_row($sorted[0]);
    }
  }
  $rpt->finish();


}
