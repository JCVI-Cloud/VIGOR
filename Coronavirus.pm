use strict;
our $myData;

{no warnings 'redefine';

# initialize VIGOR for this genome
sub customize_defaults {
	
	# parameter overrides
	set_parameter( "reference_db", "$myData/coronavirus_pep.db2" );
	set_parameter( "candidate_blastopts", "-p blastx -e 0.001 -b 2500 -F \"\"" );
	set_parameter( "min_candidate_pctsimilarity", 30 );
	set_parameter( "min_candidate_sbjcoverage", 50 );
	set_parameter( "candidate_overlap_threshold",   0.85 );
	
	set_parameter( "rawgene_reblastopts", "-p blastx -e 0.001 -v 0 -b 10 -F \"\"" );
	set_parameter( "min_rawgene_pctsimilarity", 30 );
	set_parameter( "min_rawgene_sbjcoverage", 50 );
	set_parameter( "rawgene_overlap_threshold",   0.85 );
	
	set_parameter( "min_generef_pctcoverage", 50 );
	set_parameter( "min_generef_pctsimilarity", 30 );
	set_parameter( "min_generef_pctsimilarity25", 0 );
	set_parameter( "min_generef_pctidentity", 20 );
	set_parameter( "max_generef_pcttrunc5", 50 );
	set_parameter( "max_generef_pcttrunc3", 50 );

	#  cleavage site overrides (for mature peptides)
	my $cleavage_sites;
	$$cleavage_sites{"LQ"}{offset} = 2;
	$$cleavage_sites{"LQ"}{bonus} = 10;
	$$cleavage_sites{"GG"}{offset} = 2;
	$$cleavage_sites{"GG"}{bonus} = 9;
	$$cleavage_sites{"FQ"}{offset} = 2;
	$$cleavage_sites{"FQ"}{bonus} = 9;
	$$cleavage_sites{"VQ"}{offset} = 2;
	$$cleavage_sites{"VQ"}{bonus} = 9;
	$$cleavage_sites{"MQ"}{offset} = 2;
	$$cleavage_sites{"MQ"}{bonus} = 9;
	$$cleavage_sites{"RR"}{offset} = 2;
	$$cleavage_sites{"RR"}{bonus} = 8;
	$$cleavage_sites{"QA"}{offset} = 1;
	$$cleavage_sites{"QA"}{bonus} = 10;
	$$cleavage_sites{"QR"}{offset} = 1;
	$$cleavage_sites{"QR"}{bonus} = 8;
	$$cleavage_sites{"QS"}{offset} = 1;
	$$cleavage_sites{"QS"}{bonus} = 10;
	$$cleavage_sites{"QN"}{offset} = 1;
	$$cleavage_sites{"QN"}{bonus} = 7;
	$$cleavage_sites{"GA"}{offset} = 1;
	$$cleavage_sites{"GA"}{bonus} = 5;
	$$cleavage_sites{"GR"}{offset} = 1;
	$$cleavage_sites{"GR"}{bonus} = 3;
	$$cleavage_sites{"GS"}{offset} = 1;
	$$cleavage_sites{"GS"}{bonus} = 5;
	$$cleavage_sites{"GN"}{offset} = 1;
	$$cleavage_sites{"GN"}{bonus} = 2;
	set_cleavage_sites( $cleavage_sites );
	
	return;
}

sub candidate_comparison_exception {
	my ( $genome, $ihit, $jhit ) = @_;
	
	if ( get_reference_name( $$ihit{subject_id} ) eq "PP1a"
			&& get_reference_name( $$jhit{subject_id} ) eq "PP1ab" ) {
		return 0;
	}
	elsif ( get_reference_name( $$ihit{subject_id} ) eq "PP1ab"
			&& get_reference_name( $$jhit{subject_id} ) eq "PP1a" ) {
		return 0;
	}
	return undef;
}

sub customize_gene_set {
	my ( $genome, $genes ) = @_;
	
	my $orf1ab;
	for my $gene( @$genes ) {
		if ( $$gene{is_mutation} ) {
		}
		elsif ( $$gene{gene_name} eq "PP1ab" ) {
			$orf1ab = $gene;
		}
	}
	if ( defined $$orf1ab{mature_peps}  ) {
		for my $gene( @$genes ) {
			if ( $$gene{is_mutation} ) {
			}
			elsif ( $$gene{gene_name} eq "PP1a" && defined $$gene{mature_peps} ) {
				my $orf1a = $gene;
				my @pepA = @{$$orf1a{mature_peps}};
				my @pepAB = @{$$orf1ab{mature_peps}};
				for my $a ( 0..@pepA-1 ) {
					for my $ab ( 0..@pepAB-1 ) {
						if ( $pepA[$a]{pep_start} == $pepAB[$ab]{pep_start}
								&& $pepA[$a]{pep_end} == $pepAB[$ab]{pep_end} ) {
							$pepA[$a] = undef;
							last;
						}
					}
				}
				my @tmp = remove_undefs( @pepA );
				$$orf1a{mature_peps} = \@tmp;
			}
		}
	}
}

sub before_validate_gene {
	my ( $genome, $gene ) = @_;
	
	if ( $$gene{gene_name} eq "PP1ab" ) {
		if ( @{$$gene{exons}} > 2 ) {
			my @frameshifts;
			my @exons = @{$$gene{exons}};
			for my $i ( 1..@exons-1 ) {
				my $left = $exons[$i-1]{dna_end} + $$gene{orientation};
				my $right = $exons[$i]{dna_start} - $$gene{orientation};
				if ( $left > $right ) {
					push @frameshifts, "$right..$left";
				}
				elsif ( $right > $left ) {
					push @frameshifts, "$left..$right";
				}
				else {
					push @frameshifts, $left;
				}
			}
			my $msg = "multiple frameshifts at " . join(", ", @frameshifts );
			$msg =~ s/, ([^,]*)$/ and $1/;
			if ( defined $$gene{warning} && length( $$gene{warning} ) ) {
				$$gene{warning} .= ", $msg"; 
			}
			else {
				$$gene{warning} = $msg;
			}
			$$gene{is_mutation} = 1;
			return 0;			
		}
	} 
	return undef;
}
} 
1;
