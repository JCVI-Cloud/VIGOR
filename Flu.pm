use strict;
our $myData;

{no warnings 'redefine';
# see also: VIGOR_developers_reference_guide
#
# PLEASE NOTE: code samples for the various subroutines are
#              not necessarily consistent with one another
# =======================================================================

my $vigorspace = get_parameter("vigorspace");
my $gene1 = undef;
my $rawgene1 = undef;
#===================================================================================
#
# override parameter defaults
# load reference database
sub customize_defaults {
	
	# partial genes and splicing 
#	set_parameter( "allow_partial_genes", 0 );
	#set_parameter( "allow_splicing", 2 );
	set_parameter( "spliced_gene_list", "NS2,M2" );
#	set_parameter( "small_exon_autofind", 0 );
	set_parameter( "min_exon_size", 1 );
	set_parameter( "min_intron_size", 200 );
#	set_parameter( "frameshift_detection", 1 );

	# attempt to accomodate shorter splice variants
	set_parameter( "candidate_blastopts", "-p blastx -e 1e-20 -F \"\" -v 0 -b 50" );
	set_parameter( "rawgene_reblastopts", "-p blastx -e 1e-20 -F \"\" -v 0 -b 50" );

	# candidate region identification
	set_parameter( "min_candidate_pctsimilarity", 50 );
	set_parameter( "min_candidate_sbjcoverage", 50 );
#	set_parameter( "candidate_overlap_method", 1 );
	set_parameter( "candidate_overlap_threshold", 0.85 );
	set_parameter( "candidate_overlap_ismutual", 1 );

	# rawgene identification
	set_parameter( "rawgene_extend_span",  100 );
	set_parameter( "rawgene_extend_hsp",  25 );
	set_parameter( "min_rawgene_pctsimilarity", 50 );
	set_parameter( "min_rawgene_sbjcoverage", 50 );

	# gene validation
	set_parameter( "min_generef_pctcoverage", 50 );
	set_parameter( "min_generef_pctsimilarity", 80 );
	set_parameter( "min_generef_pctidentity", 75 );
	set_parameter( "max_generef_pcttrunc5", 50 );
	set_parameter( "max_generef_pcttrunc3", 50 );
	set_parameter( "min_generef_pctsimilarity25", 85 );
	set_parameter( "detailed_exceptions", 1 );

	# use pre-built clustered reference db for candidate identification
	set_parameter( "candidate_refdb", "$myData/flu_candidates" );  	

#	cluster the reference sequences for candidate identification
#	set_parameter( "cluster_candidaterefs", 1 );  	
#	set_parameter( "candidateref_clusterL", 0.90 );		# mutual coverage	
#	set_parameter( "candidateref_clusterS", 90 );		# %identity

	# build reference database
	open( PAL, ">$vigorspace/refdb.pal" );
	print PAL "TITLE VIGOR reference sequence\n";
	print PAL "DBLIST $myData/Influenza_pep_db ";
	print PAL "$myData/PB1_F2_pep ";
	print PAL "$myData/FluB_NB_pep\n";
	close(PAL);

	set_parameter( "reference_db", "$vigorspace/refdb" );

	return;
}

#sub before_find_candidates {
#	my ( $genome ) = @_;	
#	set_parameter( "allow_splicing", 0 );
#	$gene1 = undef;
#	$rawgene1 = undef;
#}

# customized post processing of candidate regions
#sub after_find_candidates {
#	my ( $genome, $candidates ) = @_;
#	set_parameter( "allow_splicing", 1 );
#}

#sub before_find_rawgene {
#	my ( $genome, $candidate ) = @_;
#	
#	my $seq_name = get_reference_name( $$candidate{subject_id} );
##print "seqname=\"$seq_name\"\n";
#	if ( $seq_name eq "NS2" || $seq_name eq "M2" ) {
#		set_parameter( "allow_splicing", 1 );
#	}
#	else {
#		set_parameter( "allow_splicing", 0 );
#	}
#}

#sub	after_find_rawgenes {
#	my ( $genome, $candidate, $rawgenes ) = @_;
#	
#	if ( @$rawgenes > 1 ) {
#		@$rawgenes =
#			sort {
#				get_reference_name( $$a{subject_id} ) cmp get_reference_name( $$b{subject_id} )
#			} @$rawgenes;
#	}
#
#	set_parameter( "allow_splicing", 1 );
#}

# customized initialization of find_gene
# (e.g. tweak parameters for specific genomic sequence)
#sub before_find_gene {
#	my ( $genome, $rawgene ) = @_;
#
#	my $seq_name = get_reference_name( $$rawgene{subject_id} );
##print "seqname=\"$seq_name\"\n";
#	if ( $seq_name eq "NS2" || $seq_name eq "M2" ) {
#		set_parameter( "allow_splicing", 1 );
#	}
#	else {
#		set_parameter( "allow_splicing", 0 );
#		return;
#	}
#
#	if ( defined $gene1 ) {
#		if ( $$gene1{orientation} != $$rawgene{orientation} ) { return }
#		if ( $seq_name eq "NS2" && $$gene1{gene_name} eq "NS1" ) {
#		}
#		elsif ( $seq_name eq "M2" && $$gene1{gene_name} eq "M1" ) {
#		}
#		else {
#			return;
#		}
#	}
#	elsif ( defined $rawgene1 ) {
#		my $raw_name = get_reference_name( $$rawgene1{subject_id} );
#		if ( $$rawgene1{orientation} != $$rawgene{orientation} ) { return }
#		if ( $seq_name eq "NS2" && $raw_name eq "NS1" ) {
#		}
#		elsif ( $seq_name eq "M2" && $raw_name eq "M1" ) {
#		}
#		else {
#			return;
#		}
#	}
#	else {
#		return;
#	}
#			
#
#	my @hsps = sort { compare_qry_positions( $a, $b ) } @{$$rawgene{hsps}};	
##print "\nFLU-RAWGENE RAW IN\n";
##print_blasthits( $rawgene );
#	my $ori = $$rawgene{orientation};
#	my $subjleft;
#	my $dna_start;
#	if ( defined $gene1 ) {
#		my $exon = ${$$gene1{exons}}[0];
#		$subjleft = $$exon{subject_start};
#		$dna_start = $$exon{dna_start};
##print "\nFLU-RAWGENE EXON1\n";
##print_genehits( $gene1 );
#	}
#	else {
#		$subjleft = $$rawgene1{subject_left};
#		if ( $ori == 1 ) {
#			$dna_start = $$rawgene1{query_left};
#		}
#		else {
#			$dna_start = $$rawgene1{query_right};
#		}
##print "\nFLU-RAWGENE EXON1\n";
##print_blasthits( $rawgene1 );
#	}
##print "S: $subjleft  Q: $dna_start\n";
#	
## 	do we already have the first exon?
#	if ( @hsps > 1 ) {
#		if ( $$rawgene{subject_left} < 5 ) { return }
#
#		# is first exon truncated?
#		my $missing = ( $$rawgene{subject_left}-1 ) * 3 + 5;
#		if ( $ori == 1 ) {
#			if ( $$rawgene{query_left} < $missing ) { return }
#		}
#		elsif ( $$rawgene{query_right} >= length( $$genome{sequence} )-$missing ) {
#			return;
#		}
#		
#		# cannot use first exon, discard it
#		shift @{$$rawgene{hsps}};
#	}
#
#	if ( $subjleft <= 8 ) {
#		my $subjright = maxval( $$rawgene{subject_left}-1, 8 );		
#		my $dna_end = $dna_start + $ori * 3 * ( $subjright - $subjleft + 1 ) - $ori;
#		
#		# does the exon fit?
#		if ( $ori == 1 ) {
#			if ( $dna_end > $$rawgene{query_left} - 100 ) { return }
#		}
#		else {
#			if ( $dna_end < $$rawgene{query_right} + 100 ) { return }			
#		}
#		
#		my $subjcov = $subjright - $subjleft + 1;
#		
#		my $numid = maxval( 1, int( $subjcov / 5.0 + 0.5 ) );
#		my $numsim = maxval( 1, int ( $subjcov / 3.0 + 0.5 ) );
#		
#		add_new_hsp( $genome, $rawgene, $subjleft, $subjright,
#			minval( $dna_start, $dna_end ), maxval( $dna_start, $dna_end ),
#			$numid, $numsim );
#	}
##print "\nFLU-RAWGENE OUT\n";
##print_blasthits( $rawgene, @{$$rawgene{hsps}} );
##print "\n";
#
#	return;
#}

# override default selection of start codon
# return genomic position of first base of start codon (or of first peptide if partial)
#sub start_codon_exception {
#	my ( $sequence, $rawgene, $exons ) = @_;
#
#	if ( ! defined $gene1 ) { return }
#
#	my $seq_name = get_reference_name( $$rawgene{subject_id} );
#	if ( $seq_name eq "NS2" && $$gene1{gene_name} eq "NS1" ) {
#	}
#	if ( $seq_name eq "M2" && $$gene1{gene_name} eq "M1" ) {
#	}
#	else {
#		return;
#	}
#
#	return $$gene1{start_codon};
#	return undef;
#}

# customized validation of gene
#sub gene_validation_exception {
#	my ( $gene ) = @_;
#
#	if ( $$gene{error} || $$gene{is_mutation} ) {
#		return;
#	}
#	elsif ( $$gene{gene_name} eq "NS2" || $$gene{gene_name} eq "M2" ) {
#		my @exons = @{$$gene{exons}}; 
#		if ( @exons > 2 ) {
#			$$gene{is_mutation} = 1;
#			my $last = @exons[@exons-1];
#			if ( exists $$last{missing_exon} ) {
#				$$gene{warning} = "last exon is truncated";
#			}
#			else {
#				$$gene{warning} = "improper splicing";				
#			}
#		}
#	}
#}

# customized post-processing of gene
#sub after_find_gene {
#	my ( $genome, $rawgene, $gene ) = @_;
#	
#	if ( $$gene{is_mutation} ) {
#		my $seq_name = get_reference_name( $$rawgene{subject_id} );
#		if ( $seq_name eq "NS1" || $seq_name eq "M1" ) { $rawgene1 = $rawgene }
#	}
#	elsif ( $$gene{gene_name} eq "NS1" || $$gene{gene_name} eq "M1" ) {
#		$gene1 = $gene;
#	}
#
#	set_parameter( "allow_splicing", 1 );
#}

# return standardized reference name
sub get_reference_name {
	my ($reference_id) = @_;
	my $seq = get_reference_seq($reference_id);
	
	my $seq_name;
	if ( $$seq{defline} =~ /\((\w+)\)\// ) {
		$seq_name = $1;
		if ( $seq_name eq "HE" ) { $seq_name = "HA" }
	}

	if ( $$seq{defline} =~ /neuraminidase/i ) {
		$seq_name = "NA";		
	}
	elsif ( $$seq{defline} =~ /non[- ]*structural *protein *2/i ) {
		$seq_name = "NS2";		
	}
	elsif ( $$seq{defline} =~ /NSP{0,1} *2/ ) {
		$seq_name = "NS2";		
	}
	elsif ( $$seq{defline} =~ /non[- ]*structural *protein/i ) {
		$seq_name = "NS1";		
	}
	elsif ( $$seq{defline} =~ /NSP{0,1}[ 1]*/ ) {
		$seq_name = "NS1";		
	}
	
	if ( !defined $seq_name ) {
		for my $name ( "BM2", "CM", "CM1", "CM2", "HA", "HE", "MP", "M1", "M2" , "NA", "NB", "NS", "NS1", "NS2", "PA", "PB1", "PB2" ) {
			if ( " $$seq{defline} " =~ /\W$name\W/ ) {
				$seq_name = $name;
				last;
			}
		}
	}
	
	if ( $seq_name eq "NA" ) {
		if ( " $$seq{defline} " =~ /\W(CM[12]{0,1})\W/ ) {
			$seq_name = "$1";
			if ( $seq_name eq "CM" ) {
				$seq_name = "CM1";
			}
		}
	}
	elsif ( $seq_name eq "MP" ) {
		if ( " $$seq{defline} " =~ /\W(M2)\W/ ) {
			$seq_name = "M2";
		}
		elsif ( " $$seq{defline} " =~ /\W(BM2)\W/ ) {
			$seq_name = "BM2";
		}
		else {
			$seq_name = "M1";
		}
	}	
	elsif ( $seq_name eq "PB1" ) {
		if ( " $$seq{defline} " =~ /\W(F2)\W/ ) {
			$seq_name = "$seq_name-F2";
		}
	}	
	elsif ( $seq_name eq "CM" ) {
		$seq_name = "CM1";
	}	
	elsif ( $seq_name eq "NS" ) {
		$seq_name = "NS1";
	}	

	if ( ! defined $seq_name ) { $seq_name = $$seq{id} }
#print " $seq_name  |  $$seq{defline}\n";
	return $seq_name;
}

# return standardized reference definition
sub get_reference_definition {
	my ($reference_id) = @_;

	my $seq = get_reference_seq($reference_id);
	my $definition = $$seq{defline};
	$definition =~ s/\s+/ /g;
	$definition =~ s/^[^ ]* *//;
	$definition =~ s/[^ ]*\/[^ ]*\/[^ ]*\/[^ ]*\/[^ ]* //g;
	$definition =~ s/(\[.+) \(.+\)/$1/g;
	$definition =~ s/(non)[- ]*(struct)/$1-$2/gi;
	$definition =~ s/\/.+\/ //g;
	$definition =~ s/ \(.+/\]/g;

	my $seq_name = get_reference_name( $reference_id );
	my ( $main_name, $sub_name ) = split /-/, $seq_name;
	if ( ! defined $sub_name ) {
		$sub_name = "";
	}
	if ( $main_name eq "HA" ) {
		$definition =~ s/.+ \[/hemagglutinin \[/;
	}
	elsif ( $main_name eq "NA" ) {
		if ( length( $sub_name ) ) {
			$definition =~ s/.+ \[/sub_name protein \[/;
		}
		else {
			$definition =~ s/.+ \[/neuraminidase \[/;
		}
	}
	elsif ( $main_name eq "PA" ) {
		$definition =~ s/.+ \[/polymerase PA \[/;
	}
	elsif ( $main_name eq "PB1" ) {
		if ( length( $sub_name ) ) {
			$definition =~ s/.+ \[/$seq_name protein \[/;
		}
		else {
			$definition =~ s/.+ \[/polymerase PB1 \[/;
		}
	}
	elsif ( $main_name eq "PB2" ) {
		$definition =~ s/.+ \[/polymerase PB2 \[/;
	}
	elsif ( $main_name eq "M1" || $main_name eq "M2" ) {
		$definition =~ s/.+ \[/$main_name matrix protein \[/;
	}
	elsif ( $main_name =~ /^[BC]M\d$/ ) {
		$definition =~ s/.+ \[/$main_name protein \[/;
	}
	elsif ( $main_name eq "NB" ) {
		$definition =~ s/.+ \[/NB protein \[/;
	}
	elsif ( $main_name eq "NS2" ) {
		$definition =~ s/.+ \[/non\-structural protein 2 \[/;
	}
	elsif ( $main_name eq "NS1" ) {
		$definition =~ s/.+ \[/non\-structural protein 1 \[/;
	}
	elsif ( $main_name eq "NP" ) {
		$definition =~ s/.+ \[/nucleoprotein \[/;
	}

	if ( $$seq{defline} =~ /(H\d+)N\d+/ ) {
		$definition .= " $1";
	}

	return $definition;
}
}
1;