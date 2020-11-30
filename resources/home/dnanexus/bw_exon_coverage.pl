#!/bin/env perl
# report average coverage for coding exon bases from BigWig for a
# set of genes/transcripts
#
# Michael Edmonson 11/2013 -

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;

use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;

use FileUtils qw(read_simple_file);
use Reporter;
use MiscUtils qw(median average dump_die);
use Chr2GeneParser;
use TdtConfig;
use SJPreferredIsoform;

my %FLAGS;
my @BW_FILES;
my @EXON_FILTER;
my @EXTRA_GENES;

my $CHR2GENE_FEATURE = "CodingExon";

my $CLINICAL_MIN_COVERAGE = 20;
my $CLINICAL_MIN_UNCOVERED_EXON_PERCENT_TO_REPORT = 5;

GetOptions(\%FLAGS,
	   "-genome=s",
	   "-gene-list=s",
	   # list of gene symbols to restrict to
	   # (note that if used there may well be multiple isoforms)

	   "-isoform-list=s",
	   # list of NM_ accessions to restrict to

	   "-bw=s" => \@BW_FILES,
	   "-bw-list=s",

	   "-chr2gene=s",
	   # e.g. /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/chr2gene/

	   "-out=s",
	   # manual outfile (if only one input .bw)

	   "-suffix=s",
	   # suffix to apply to outfile

	   "-combine",

	   "-exon=i" => \@EXON_FILTER,

	   "-cache",

	   "-germline-reviewable",
	   "-germline-reportable",

	   "-debug-parse=s",

#	   "-feature=s" => \$CHR2GENE_FEATURE,
	   # change to Exon to include non-coding exons as well
	   # FIX ME: should this be a list, i.e. both Exon and CodingExon?
	   # not sure

	   "-clinical-summary",
	   # brief report of low-coverage exons
	   "-clinical-min-coverage=i" => \$CLINICAL_MIN_COVERAGE,
	   "-clinical-min-report-percent=f" => \$CLINICAL_MIN_UNCOVERED_EXON_PERCENT_TO_REPORT,

	   "-preferred-only",

	   "-generate-clinical-tests",

	   "-dump-exon-counts",
	   "-root=s",
	   "-thresholds=s",

	   "-add-gene=s" => \@EXTRA_GENES,
	  );

if ($FLAGS{"generate-clinical-tests"}) {
  generate_clinical_tests();
  exit(0);
}

my %QC_IGNORE_GENES = map {$_, 1} qw(
				      TERC
				   );
# don't complain about missing non-coding genes

my %WANTED_EXONS = map {$_, 1} @EXON_FILTER;

if (my $fn = $FLAGS{"bw-list"}) {
  my $list = read_simple_file($fn);
  push @BW_FILES, @{$list};
}

die "specify -bw or -bw-list [file]" unless @BW_FILES;

my ($gene_list, $iso_list);

if (my $f = $FLAGS{"gene-list"}) {
  # restrict isoforms to just those for specified gene list
  $gene_list = read_simple_file($f);
} elsif (my $f2 = $FLAGS{"isoform-list"}) {
  # restrict isoforms to specific accessions
  $iso_list = read_simple_file($f2);
} elsif ($FLAGS{"germline-reviewable"}) {
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $file = $config_genome->{CLINCLS_CANCER_RELATED_GENES_FILE} || die;
  open(TMP, $file) || die;
  my %all_reviewable;
  while (<TMP>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 2;
    $all_reviewable{$f[0]} = 1;
  }
  $gene_list = [sort keys %all_reviewable];
  printf STDERR "debug: wanted genes: %s\n", join ", ", @{$gene_list};
  die "germline reviewable still needs work";
} elsif ($FLAGS{"germline-reportable"}) {
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
#  my $file = $config_genome->{CLINCLS_GERMLINE_REPORTABLE_GENES} || die;
  my $file = $config_genome->{CLINCLS_GERMLINE_REPORTABLE_GENES_CHR2GENE} || die;

  my $list = read_simple_file($file);
  push @{$list}, @EXTRA_GENES if @EXTRA_GENES;
  $gene_list = [ sort @{$list} ];
  printf STDERR "debug: wanted genes: %s\n", join ", ", @{$gene_list};
} else {
  #  die "specify -gene-list [file] or -isoform-list [file]\n";
  # global report: still needs work, may fail due to multiple chr mappings etc.
}

my $isoforms = load_isoforms(
			     "-genes" => $gene_list,
			     "-isoforms" => $iso_list,
			    );

my @headers = (
	       "chr",
	       "gene",
	       "strand",
	       "accession",
	       "# bases >=20x",
	       "# bases <20x",
	       '%covered >=20x',
	       '%covered >=30x',
	       '%covered >=40x',
	      );

if (@EXON_FILTER) {
  push @headers, (
		  map {"Exon" . $_} @EXON_FILTER
		 );
} else {
  push @headers, (
		  map {"Exon" . $_} (1..125)
		 );
}

if ($FLAGS{combine}) {
  generate_combined_report();
  exit(0);
} elsif ($FLAGS{"clinical-summary"}) {
  generate_clinical_report();
  # a variant of combined report
  exit(0);
}

my %saw;
foreach my $fn (@BW_FILES) {
  # unique outfile check:
  # do this BEFORE processing so we don't imply processing was successful
  # by producing some output before crashing
  next unless $fn =~ /\w/;
  my $bn = basename($fn);
  die "ERROR: non-unique basename for $fn!" if $saw{$bn};
  $saw{$bn} = 1;
}

foreach my $fn (@BW_FILES) {
  next unless $fn =~ /\w/;
  # in case blank lines in file list
  die "where is $fn?" unless -s $fn;

  my $outfile;
  if ($outfile = $FLAGS{out}) {
    die "can't specify -out with multiple -bw" if @BW_FILES > 1;
  } elsif (my $s = $FLAGS{suffix}) {
    $outfile = sprintf '%s.%s.coverage.tab', basename($fn), $s;
  } else {
    $outfile = sprintf '%s.coverage.tab', basename($fn);
  }
  next if $FLAGS{cache} and -s $outfile;

  printf STDERR "opening %s...\n", $fn;
  my $bw = Bio::DB::BigFile->bigWigFileOpen($fn) || die;

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => \@headers,
			);

  foreach my $accession (sort keys %{$isoforms}) {
    print STDERR "processing $accession...\n";
    my $exons = $isoforms->{$accession};

    my $chrom = get_singleton($exons, "chrom");
    my $gene = get_singleton($exons, "gene");
    my $strand = get_singleton($exons, "strand");

    my %r;
    $r{chr} = $chrom;
    $r{gene} = $gene;
    $r{strand} = $strand;
    $r{accession} = $accession;

    my $total_bases = 0;
    my $total_bases_less_20 = 0;
    my $total_bases_ge_20 = 0;
    my $total_bases_ge_30 = 0;
    my $total_bases_ge_40 = 0;

    foreach my $exon (@{$exons}) {
      my $start = $exon->{start_inbase} || die;
      my $end = $exon->{end} || die;
      # 1-based / in-base

      my $exon_number = $exon->{exon_number} || die;

      next if @EXON_FILTER and not($WANTED_EXONS{$exon_number});

      my $intervals = $bw->bigWigIntervalQuery($chrom, $start - 1 => $end);
      my @values;
      for (my $i=$intervals->head;$i;$i=$i->next) {
	my $start = $i->start;
	my $end   = $i->end;
	my $val   = $i->value;
	#	    printf "%s\n", join " ", $start, $end, $val;
	$total_bases++;
	push @values, $val;
	if ($val < 20) {
	  printf STDERR "only %d reads at %s %s.%d\n", $val, $fn, $chrom, $end;
	  $total_bases_less_20++;
	} else {
	  $total_bases_ge_20++;
	  $total_bases_ge_30++ if $val >= 30;
	  $total_bases_ge_40++ if $val >= 40;
	}
      }
      printf STDERR "%s %s exon %d %d-%d total:%d\n", $gene, $accession, $exon_number, $start, $end, $total_bases;

      #	$r{"Exon" . $exon_number} = median(\@values);
      my $key = "Exon" . $exon_number;
      die "say what now?" if exists $r{$key};
      $r{$key} = int(average(\@values));
    }
    $r{"# bases <20x"} = $total_bases_less_20;
    $r{"# bases >=20x"} = $total_bases_ge_20;

    $r{'%covered >=20x'} = sprintf '%d%%', $total_bases_ge_20 * 100 / $total_bases;
    $r{'%covered >=30x'} = sprintf '%d%%', $total_bases_ge_30 * 100 / $total_bases;
    $r{'%covered >=40x'} = sprintf '%d%%', $total_bases_ge_40 * 100 / $total_bases;

    foreach (@headers) {
      $r{$_} = "" unless defined $r{$_};
    }
    $rpt->end_row(\%r);
  }

  $rpt->finish();
}

sub load_isoforms {
  my (%options) = @_;

  my $dir = $FLAGS{chr2gene};
  unless ($dir) {
    my $genome = $FLAGS{genome} || die "-genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $dir = $config_genome->{CHR2GENE} || die;
  }
  die "no chr2gene dir" unless $dir;
  printf STDERR "chr2gene: %s\n", $dir;

#  my @files = glob($dir . "/*.txt");
  my @files = glob($dir . "/chr*.txt");
  my $p = new Chr2GeneParser();

  if (my $l = $options{"-genes"}) {
    $p->restrict_genes($l);
  } elsif (my $l2 = $options{"-isoforms"}) {
    $p->restrict_isoforms($l2);
  }
  $p->restrict_feature_types([ $CHR2GENE_FEATURE ]);
  $p->bucket_by_isoform(1);

  my %isoforms;

  my $debug_parse = $FLAGS{"debug-parse"};
  foreach my $fn (@files) {
    next if -l $fn;
    # e.g. "chr23.gene.txt" is symlinked to chrX.gene.txt
    next if $debug_parse and $fn !~ /$debug_parse/;
    print STDERR "$fn...\n";
    $p->parse("-file" => $fn);
  }
  $p->qc_wanted_genes("-ignore" => \%QC_IGNORE_GENES);

  my $all_isoforms = $p->all_isoforms();

  if ($FLAGS{"preferred-only"}) {
    #
    # filter isoforms to preferred isoform for each gene only
    #
    my $genome = $FLAGS{genome} || die "-genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    my $matrix = $config_genome->{GENE_TRANSCRIPT_MATRIX} || die;
    my $sjpi = new SJPreferredIsoform("-file" => $matrix);

    my %gene2isoform;
    foreach my $i (keys %{$all_isoforms}) {
      my %g;
      foreach my $thing (@{$all_isoforms->{$i}}) {
	my $gene = $thing->{gene} || die;
	$g{$gene} = 1;
      }
      die unless scalar keys %g == 1;
      foreach my $g (keys %g) {
	push @{$gene2isoform{$g}}, $i;
      }
    }

    my @toss;
    foreach my $gene (keys %gene2isoform) {
      my @all_i = @{$gene2isoform{$gene}};
      my $pref = $sjpi->get_preferred_isoform($gene);
      if ($pref) {
	my $found_pref;
	foreach my $i (@all_i) {
	  if ($i eq $pref) {
	    $found_pref = $i;
	  } else {
	    push @toss, $i;
	  }
	}
      } else {
	# might not have one due to gene symbol incompatibility between
	# (A) reportable gene list
	# (B) symbol as found in chr2gene
	# (C) symbol as found in preferred isoforms file
	die "no preferred isoform for $gene " . scalar @all_i if @all_i > 1;
	# only a problem if > 1 isoform for the gene in question
      }

    }
    delete @{$all_isoforms}{@toss};
  }

  if ($FLAGS{"dump-exon-counts"}) {
    my $exons = 0;
    foreach my $i (keys %{$all_isoforms}) {
      foreach my $ref (@{$all_isoforms->{$i}}) {
	$exons++;
      }
    }
    printf STDERR "isoforms:%s exons:%d\n",
      scalar(keys %{$all_isoforms}),
	$exons;
    exit(0);
  }

  return $all_isoforms;
}

sub get_singleton {
  # "should" be only one value (make sure)
  my ($list, $key) = @_;
  my %things = map {$_->{$key}, 1} @{$list};
  die sprintf "ERROR: multiple hits for %s (%s)", $key, join(",", keys %things) if keys %things > 1;
  my ($thing) = keys %things;
  return $thing;
}

sub generate_combined_report {
  # 12/10/2013:
  # new report combining coverage values from multiple files.
  my $outfile = $FLAGS{out};

  my %file2bw;
  foreach my $fn (@BW_FILES) {
    next unless $fn =~ /\w/;
    # in case blank lines in file list
    die "where is $fn?" unless -s $fn;
    $file2bw{$fn} = Bio::DB::BigFile->bigWigFileOpen($fn) || die;
    unless ($outfile) {
      my @things = basename($fn);
      push @things, $FLAGS{suffix} if $FLAGS{suffix};
      push @things, "coverage", "tab";
      $outfile = join ".", @things;
    }
  }
  my $total_bw = scalar keys %file2bw;

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => \@headers,
			);

  foreach my $accession (sort keys %{$isoforms}) {
    print STDERR "processing $accession...\n";
    my $exons = $isoforms->{$accession};

    my $chrom = get_singleton($exons, "chrom");
    my $gene = get_singleton($exons, "gene");
    my $strand = get_singleton($exons, "strand");

    my %r;
    $r{chr} = $chrom;
    $r{gene} = $gene;
    $r{strand} = $strand;
    $r{accession} = $accession;

    my $total_bases = 0;
    my $total_bases_less_20 = 0;
    my $total_bases_ge_20 = 0;
    my $total_bases_ge_30 = 0;
    my $total_bases_ge_40 = 0;

    my %saw_exons;

    foreach my $exon (@{$exons}) {
      my $ex_start = $exon->{start_inbase} || die;
      my $ex_end = $exon->{end} || die;
      # 1-based / in-base
      my $exon_number = $exon->{exon_number} || die;
      die "duplicate entry for $exon_number" if $saw_exons{$exon_number};
      $saw_exons{$exon_number}++;

      my @values = ();

      my $bw_num = 0;
      foreach my $fn (keys %file2bw) {
	# accumulate total coverage across all .bw files
	my $bw = $file2bw{$fn};
	$bw_num++;
	die "WTF" if $bw_num > $total_bw;

	my $intervals = $bw->bigWigIntervalQuery($chrom,
						 $ex_start - 1,
						 $ex_end);
	my $ai = 0;
	for (my $i=$intervals->head; $i; $i=$i->next) {
	  my $start = $i->start;
	  my $end   = $i->end;
	  $values[$ai] += $i->value;

#	  printf STDERR "debug exon %d range=%d-%d, %s.%s in %s coverage=%d\n", $exon_number, $ex_start, $ex_end, $chrom, $end, $fn, $values[$ai];

	  if ($values[$ai] < 20 and $bw_num == $total_bw) {
	    printf STDERR "only %d total reads at %s.%d in %s\n", $values[$ai], $chrom, $end, join ",", sort keys %file2bw;
	  }

	  $ai++;
	}
      }

      foreach my $val (@values) {
	$total_bases++;
	if ($val < 20) {
	  $total_bases_less_20++;
	} else {
	  $total_bases_ge_20++;
	  $total_bases_ge_30++ if $val >= 30;
	  $total_bases_ge_40++ if $val >= 40;
	}
      }

      printf STDERR "%s %s exon %d %d-%d total:%d\n", $gene, $accession, $exon_number, $ex_start, $ex_end, $total_bases;

      #	$r{"Exon" . $exon_number} = median(\@values);
      my $key = "Exon" . $exon_number;
      die "say what now?" if exists $r{$key};
      $r{$key} = int(average(\@values));
    }
    $r{"# bases <20x"} = $total_bases_less_20;
    $r{"# bases >=20x"} = $total_bases_ge_20;

    $r{'%covered >=20x'} = sprintf '%d%%', $total_bases_ge_20 * 100 / $total_bases;
    $r{'%covered >=30x'} = sprintf '%d%%', $total_bases_ge_30 * 100 / $total_bases;
    $r{'%covered >=40x'} = sprintf '%d%%', $total_bases_ge_40 * 100 / $total_bases;

    foreach (@headers) {
      $r{$_} = "" unless defined $r{$_};
    }
    $rpt->end_row(\%r);
  }

  $rpt->finish();
}

sub generate_clinical_report {
  # 9/3/2015:
  # brief summary of low-coverage exons for combined WGS/WES coverage.
  my $outfile = $FLAGS{out};

  die "multiple files expected (WGS/WES)" unless @BW_FILES >= 2;

  my %file2bw;
  foreach my $fn (@BW_FILES) {
    next unless $fn =~ /\w/;
    # in case blank lines in file list
    die "where is $fn?" unless -s $fn;
    $file2bw{$fn} = Bio::DB::BigFile->bigWigFileOpen($fn) || die;
    unless ($outfile) {
      my @things = basename($fn);
      push @things, sprintf "mincov_%d", $CLINICAL_MIN_COVERAGE;
      push @things, sprintf "minpct_%02d", $CLINICAL_MIN_UNCOVERED_EXON_PERCENT_TO_REPORT;
      $outfile = sprintf "%s.tab", join "_", @things;
    }
  }
  my $total_bw = scalar keys %file2bw;

  my @headers = qw(
		    min_coverage
		    min_low_coverage_fraction_to_report
		    gene
		    accession
		    strand
		    exon
		    chr
		    start
		    end
		    median_coverage
		    low_coverage_fraction
		 );

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => \@headers,
			);
  $rpt->write_headers();

  foreach my $accession (sort keys %{$isoforms}) {
    print STDERR "processing $accession...\n";
    my $exons = $isoforms->{$accession};

    my $chrom = get_singleton($exons, "chrom");
    my $gene = get_singleton($exons, "gene");
    my $strand = get_singleton($exons, "strand");

    my %saw_exons;

    foreach my $exon (@{$exons}) {
      my $ex_start = $exon->{start_inbase} || die;
      my $ex_end = $exon->{end} || die;
      # 1-based / in-base
      my $exon_number = $exon->{exon_number} || die;
      die "duplicate entry for $exon_number" if $saw_exons{$exon_number};
      $saw_exons{$exon_number}++;

      my @values = ();

      my $bw_num = 0;
      #
      # accumulate total coverage across all .bw files:
      #
      foreach my $fn (keys %file2bw) {
	my $bw = $file2bw{$fn};
	$bw_num++;
	die "WTF" if $bw_num > $total_bw;

	my $intervals = $bw->bigWigIntervalQuery($chrom,
						 $ex_start - 1,
						 $ex_end);
	my $ai = 0;
	for (my $i=$intervals->head; $i; $i=$i->next) {
	  my $start = $i->start;
	  my $end   = $i->end;
	  $values[$ai] += $i->value;

#	  printf STDERR "debug exon %d range=%d-%d, %s.%s in %s coverage=%d\n", $exon_number, $ex_start, $ex_end, $chrom, $end, $fn, $values[$ai];

	  if ($values[$ai] < $CLINICAL_MIN_COVERAGE and $bw_num == $total_bw) {
	    printf STDERR "only %d total reads in %s at %s.%d in %s\n",
	      $values[$ai], $gene, $chrom, $end, join ",", sort keys %file2bw;
	    # debug for manual QC
	  }

	  $ai++;
	}
      }

      my $total_bases = 0;
      my $low_bases = 0;
      foreach my $val (@values) {
	$total_bases++;
	$low_bases++ if $val < $CLINICAL_MIN_COVERAGE;
      }
      my $fraction_low = $low_bases / $total_bases;

      my $cutoff = $CLINICAL_MIN_UNCOVERED_EXON_PERCENT_TO_REPORT / 100;

      my $fail = $fraction_low >= $cutoff;

      printf STDERR "%s %s exon %d %s, loc:%s:%d-%d median=%d, low_fraction=%.3f cutoff=%.2f\n",
	$gene, $accession, $exon_number,
	  ($fail ? "fail" : "OK"),
	    $chrom, $ex_start, $ex_end,
	      median(\@values), $fraction_low, $cutoff;

      if ($fail) {
	my %r;
	$r{min_coverage} = $CLINICAL_MIN_COVERAGE;
	$r{min_low_coverage_fraction_to_report} = $cutoff;
	$r{low_coverage_fraction} = sprintf '%.3f', $fraction_low;
	$r{gene} = $gene;
	$r{accession} = $accession;
	$r{strand} = $strand;
	$r{chr} = $chrom;
	$r{exon} = $exon_number;
	$r{start} = $ex_start;
	$r{end} = $ex_end;
	$r{median_coverage} = median(\@values);
	$rpt->end_row(\%r);
      }
    }
  }

  $rpt->finish();
}

sub generate_clinical_tests {
  my $root = $FLAGS{root} || "/cgs01/clingen/dev/tartan/index/data/ClinicalPilot/ClinicalPilot/";

  my @samples = glob("$root/SJ*");
  die unless @samples;

  my @thresholds;
  if (my $t = $FLAGS{thresholds}) {
    @thresholds = split /,/, $t;
  } else {
    @thresholds = qw(
		    5
		    10
		    20
		    30
		    40
		    50
		   );
  }

  my @cmds;
  foreach my $dir (@samples) {
    my $bn = basename($dir);
#    next unless $bn =~ /_[A-Z]$/;

    my ($bw_wgs) = glob($dir . "/WHOLE_GENOME/coverage/*.bw");
    my ($bw_exon) = glob($dir . "/EXOME/coverage/*.bw");
    if ($bw_exon and $bw_wgs) {
      foreach my $percent (@thresholds) {
	my $cmd = sprintf 'bw_exon_coverage.pl -preferred-only -bw %s -bw %s -genome GRCh37-lite -germline-reportable -clinical-summary -clinical-min-report-percent %d',
	$bw_wgs,
	  $bw_exon,
	    $percent;
	push @cmds, $cmd;
      }
    }
  }

  foreach (@cmds) {
    print "$_\n";
  }

}
