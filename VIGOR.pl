#!/usr/local/bin/perl -w
use strict;
use Getopt::Std;
use Cwd 'realpath';
use File::Basename;

my $program = realpath($0);
our $myBin = dirname($program);
our $myData = "$myBin/data";
our $myConf = "$myBin/conf";

require "$myBin/VIGOR.pm";

*STDERR = *STDOUT;
$|++;

print "initializing VIGOR...\n";
set_parameter( "command_line", $0 . " " . join( " ", @ARGV ) );
initialize_defaults();

# get command line parameters
#our ( $opt_a, $opt_c, $opt_f, $opt_G, $opt_i, $opt_l, $opt_O, $opt_p, $opt_s, $opt_t, $opt_v, $opt_x );
#getopt('acfGilOpstvx');
my %args;
&getopts( 'a:c:f:G:i:l:O:p:s:t:v:x:', \%args );
my ( $opt_a, $opt_c, $opt_f, $opt_G, $opt_i, $opt_l, $opt_O, $opt_p, $opt_s, $opt_t, $opt_v, $opt_x ) =
	( $args{a}, $args{c}, $args{f}, $args{G}, $args{i}, $args{l}, $args{O}, $args{p}, $args{s}, $args{t}, $args{v}, $args{x} );
if ( ! defined $opt_v ) { $opt_v = 0 }
my $verbose = $opt_v;
set_parameter( "verbose", $verbose );

if ( defined $opt_G ) {
	if ( $opt_G !~ /\.gbk/i ) { die "-G $opt_G must specify a genbank file.\n" }
	set_parameter( "reference_db", $opt_G );
	set_parameter( "max_aa_gap", 15 );
	set_parameter( "allow_splicing", 2 );
	set_parameter( "min_intron_size", 50 );
	set_parameter( "small_exon_autofind", 1 );
	set_parameter( "splicesite_search_range", 60 );
	set_parameter( "candidate_blastopts", "-p blastx -G 8 -E 2 -M BLOSUM80 -e 1e-10 -b 100 -F \"\" -g F" );
	set_parameter( "candidate_overlap_method", 1 );
	set_parameter( "candidate_overlap_ismutual", 1 );
	set_parameter( "candidate_overlap_threshold", 1.0 );
	set_parameter( "rawgene_extend_span", 0 );
	set_parameter( "rawgene_extend_hsp", 0 );
	set_parameter( "rawgene_reblastopts", "-p blastx -G 8 -E 2 -M BLOSUM80 -e 1e-10 -v 0 -F \"\" -g F" );
	set_parameter( "rawgene_overlap_threshold", 1.0 );
	set_parameter( "startcodon_search_aarange", "50%:50%" );
	set_parameter( "stopcodon_search_aaextension", "50%" );
	set_parameter( "min_generef_pctsimilarity25", 0 );
}

# load customization module
if ( defined $opt_l ) {
	require "$opt_l";
	customize_defaults();
}

# load customization configuration
if ( defined $opt_x ) {
	read_config_file( $opt_x );
}
show_parameters(); 

# parse (optional) lists of exception genes (splicing/ribosomal slippage/stop readthrus)
parse_exception_lists();

# load reference sequences
load_reference_db();

# open input and output files
open( INPUT, "<$opt_i" ) || die "cannot open the file -i \"$opt_i\"\.$!\n";
my $INPUT = *INPUT;

my $ALIGN;
if ( ! defined $opt_a && defined $opt_O ) {
	$opt_a = "$opt_O.aln";
}
if ( defined $opt_a ) {
	open( ALIGN, ">$opt_a" ) || die "cannot write to the file -a \"$opt_a\"\.$!\n";
	$ALIGN = *ALIGN;
}
if ( ! defined $ALIGN ) { print "WARNING: no path specified for alignment file (-a), file will not be generated\n" }

my $FS;
if ( ! defined $opt_f && defined $opt_O ) {
	$opt_f = "$opt_O.fs";
}
if ( defined $opt_f ) {
	open( FS, ">$opt_f" ) || die "cannot write to the file -f \"$opt_f\"\.$!\n";
	$FS = *FS;
}
if ( ! defined $FS ) { print "WARNING: no path specified for frameshift file (-f), file will not be generated\n" }

my $TBL;
if ( ! defined $opt_t && defined $opt_O ) {
	$opt_t = "$opt_O.tbl";
}
if ( defined $opt_t ) {
	open( TBL, ">$opt_t" ) || die "cannot write to the file -t \"$opt_t\"\.$!\n";
	$TBL = *TBL;
}
if ( ! defined $TBL ) { print "WARNING: no path specified for TBL file (-t), file will not be generated\n" }

my $CDS;
if ( ! defined $opt_c && defined $opt_O ) {
	$opt_c = "$opt_O.cds";
}
if ( defined $opt_c ) {
	open( CDS, ">$opt_c" ) || die "cannot write to the file -c \"$opt_c\"\.$!\n";
	$CDS = *CDS;
}
if ( ! defined $CDS ) { die "no path specified for CDS file (-c)\n" }

my $PEP;
if ( ! defined $opt_p && defined $opt_O ) {
	$opt_p = "$opt_O.pep";
}
if ( defined $opt_p ) {
	open( PEP, ">$opt_p" ) || die "cannot write to the file -p \"$opt_p\"\.$!\n";
	$PEP = *PEP;
}
if ( ! defined $PEP ) { die "no path specified for peptide file (-p)\n" }

my $STATS;
if ( ! defined $opt_s && defined $opt_O ) {
	$opt_s = "$opt_O.stats";
}
if ( defined $opt_s ) {
	open( STATS, ">$opt_s" ) || die "cannot write to the file -s \"$opt_s\"\.$!\n";
	$STATS = *STATS;
	stats_header( $STATS, undef );
}
if ( ! defined $STATS ) { print "WARNING: no path specified for statistics file (-s), file will not be generated\n" }

# for each input sequence
while ( my $genome = next_sequence( $INPUT ) ) {
	my $gene_num = 0;
	my $used_genes;
	my $fragcount = 0;
	
	while ( my $frag = next_frag( $genome ) ) {
	$fragcount++;
	my $frag_begin = $$frag{frag_offset}+1;
	my $frag_end = $$frag{frag_offset}+length($$frag{sequence});
	print "\ngenomic sequence $$genome{id}, length=". length( $$genome{sequence} ) . "  fragment $frag_begin-$frag_end";
	if ( $verbose ) { print "\n" }

	# find candidate coding regions
	before_find_candidates( $frag );
	my @candidates = find_candidates( $frag );
	after_find_candidates( $frag, \@candidates );
	print "  candidate regions: " . @candidates . "\n";

	# find raw genes in candidate regions
	for my $candidate ( @candidates ) {
		before_find_rawgenes( $frag, $candidate );
		my @rawgenes = find_rawgenes( $frag, $candidate );
		after_find_rawgenes( $frag, $candidate, \@rawgenes );
		my $region_name = get_reference_name( $$candidate{subject_id} );
		print "    candidate region: $region_name  location: " . get_blasthit_location( $candidate ) . "  potential genes: " . @rawgenes . "\n";
		# find gene (or mutation) in raw gene
		for my $rawgene ( @rawgenes ) {
			
			my @permutations = ( $rawgene );
			if ( exists $$rawgene{permutations} ) { push @permutations, @{$$rawgene{permutations}} }
			my $rawgene_name = get_reference_name( $$rawgene{subject_id} );
			
			screen_rawgene_permutations( $frag, \@permutations );
			if ( $verbose ) {
				print "\nPERMUTATIONSS\n";
				print_blasthits( @permutations );
			}

			print "      potential gene: $rawgene_name  location: " . get_blasthit_location( $rawgene ) . "  permutations: " . @permutations;
			if ( $verbose || @permutations < 2 ) { print "\n" }
			
			$gene_num++;
			my @selected;
			my $valid_count = 0;
			my $topquality;

			my $permuteid = 0;
			for my $permutation ( @permutations ) {
				$permuteid++;
				if ( ! $verbose && @permutations > 1 ) { print ".." }
				if ( $verbose ) {
					print "PERMUTATION\n";
					print_blasthits( $permutation, @{$$permutation{hsps}} );
				}
				$$permutation{gene_num} = $gene_num;
				before_find_gene( $frag, $permutation );
				my $tmpgene = find_gene( $frag, $permutation );
				after_find_gene( $frag, $permutation, $tmpgene );
				
				validate_gene( $frag, $permutation, $tmpgene );
				if ( $verbose ) {
					print "TMPGENE\n";
					print_genehits( $tmpgene );
				}

				if ( ! $verbose && @permutations > 1 ) { print "#$permuteid" }

				if ( ! defined $topquality || $$tmpgene{gene_quality} > $topquality ) {
					@selected = ( $tmpgene );
					$topquality = $$tmpgene{gene_quality};
				}
				elsif ( $$tmpgene{gene_quality} == $topquality && ! $$tmpgene{is_mutation} ) {
					push @selected, $tmpgene
				}
				if ( $verbose && @permutations > 1 ) {
					print "\nMATCH QUALITY = $$tmpgene{match_quality}  SPLICING QUALITY = $$tmpgene{splicing_quality}  GENE QUALITY = $$tmpgene{gene_quality}\n";
					print_genehits( $tmpgene );
				}
				if ( ! $$tmpgene{is_mutation} ) { $valid_count++ }

			}
			if ( ! $verbose && @permutations > 1 ) { print "\n" }
			
			# remove duplicated genes
			my @genes;
			for my $gene ( sort { $$b{gene_quality} <=> $$a{gene_quality} } @selected ) {
				if ( $valid_count > 0 && $$gene{is_mutation} ) { next }
				my $locus = get_gene_location( $gene );
				if ( $$gene{stop_truncated} ) {
					if ( $$gene{gene_name} ne $$gene{ref_id} ) { $locus .= "|$$gene{gene_name}" }
				}
				if ( defined $$used_genes{$locus} ) {
					if ( $$gene{gene_quality} > $$used_genes{$locus}{gene_quality} ) {
						for my $attribute ( gene_ref_attributes() ) {
							$$used_genes{$locus}{$attribute} = $$gene{$attribute};
						}
						if ( $verbose ) {
							print "\ngene updates existing gene\n";
							print_genehits( $$used_genes{$locus} );
						}
						next;
					}
					else {
						if ( $verbose ) {
							print "\nduplicate gene\n";
						}
						next;
					}
				}
				$$used_genes{$locus} = $gene;
				push @genes, $gene;
				if ( $$gene{is_mutation} ) { last }
			}
			
			# check for alternate splicings
			if ( @genes == 1 ) {
				print "        found " . @genes . " gene\n";
			}
			else {
				print "        found " . @genes . " genes\n";
			}
			if ( ! @genes ) {
				$gene_num--;
				next;
			}
			elsif ( @genes > 1 ) {
				for my $i ( 0..@genes-1 ) {
					$genes[$i]{gene_id} .= substr( "abcdefghijlkmnopqrstuvwxyz", $i, 1 );
					if ( length( $genes[$i]{warning} ) ) {
						$genes[$i]{warning} .= ", multiple splicings yielded equivalent protein";
					}
					else {
						$genes[$i]{warning} = "multiple splicings yielded equivalent protein";
					}
				}
			}

			# post-processing of genes
			for my $gene ( @genes ) {
				# adjust coordinates for fragment offsets
				adjust_fragment_coordinates( $frag, $gene );
	
				if ( $verbose ) {
					print "\nGENE $$gene{gene_id}\n";
					print_genehits( $gene);
				}
			}
		
			if ( $verbose && @permutations > 1 ) {
				print "\nSELECTED PERMUTATIONS\n";
				print_genehits( @genes );
			}
		}
	}
	}
	my @genes = sort { compare_gene_ids( $a, $b ) } values %$used_genes;
	if ( $fragcount > 1 ) {
		print "\nmerging gene fragments...\n";
		merge_partial_genes( $genome, \@genes );
	}

	# annotate mature peptides
	my $printmsg = 1;
	for my $gene ( @genes ) {
		if ( ! $$gene{is_mutation} ) {
			my $ref = get_reference_seq( $$gene{ref_id} );
			if ( $$ref{length} > get_parameter( "polyprotein_threshold" ) ) {
				if ( $printmsg ) {
					print "\nannotating mature peptides...\n";
					$printmsg = 0;
				}
				$$gene{mature_peps} = mapto_polyprotein($gene);
			}
			else {
				$$gene{mature_peps} = get_mature_peptides($gene);
				if ( defined $$gene{mature_peps} && $printmsg ) {
					print "\nannotating mature peptides...\n";
					$printmsg = 0;
				}
			}
		}
	}
	
	# allow for genome-specific customization of genes
	# (e.g. use orf1a gene to convert orf1b gene to orf1ab)	
	customize_gene_set( $genome, \@genes );

	print "\n";
	stats_header( undef, *STDOUT );
	stats_report( $STATS, $genome, \@genes );
	print "\n";
	print_genehits( @genes );
	print "\n";

	# we have the final set of genes for this genome 	
	# report genes
	before_report_gene( $genome, \@genes, $TBL, $CDS, $PEP, $ALIGN, $FS );
	my ( $stdt, $stdc, $stdp, $stda ) = report_gene_exception( $TBL, $CDS, $PEP, $ALIGN, $FS, $genome, \@genes );
	if ( $stdt ) { std_tbl_report( $TBL, $genome, \@genes ) }
	if ( $stdc ) { std_cds_report( $CDS, $genome, \@genes ) }
	if ( $stdp ) { std_pep_report( $PEP, $genome, \@genes ) }
	if ( $stdp ) { std_align_report( $ALIGN, $FS, $genome, \@genes ) }
	after_report_gene( $TBL, $CDS, $PEP, $ALIGN, $FS, $genome, \@genes );
}

# completed
# close files and clean up workspace
if ( defined $INPUT ) { close $INPUT }
if ( defined $STATS ) { close $STATS }
if ( defined $TBL ) { close $TBL }
if ( defined $CDS ) { close $CDS }
if ( defined $PEP ) { close $PEP }
if ( defined $ALIGN ) { close $ALIGN }
if ( defined $FS ) { close $FS }

my $retry = 0;
my $vigorspace = get_parameter( "vigorspace" );
#while ( ( system "rm -rf $vigorspace" ) && $retry++ < 10 ) {
#	sleep 10;
#}
exit(0);

