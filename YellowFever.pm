use strict;
{no warnings 'redefine';

sub customize_defaults {

	# default parameter overrides
	set_parameter( "reference_db", "/usr/local/projects/GCV/Shiliang.dir/YellowFeverV.dir/YFV_polyprotein.fasta");
	set_parameter( "min_generef_pctcoverage", 95 );
	set_parameter( "min_generef_pctsimilarity", 50 );
	set_parameter( "min_generef_pctsimilarity25",  65 );
	set_parameter( "min_generef_pctidentity", 0 );
	set_parameter( "max_generef_pcttrunc5", 5 );
	set_parameter( "max_generef_pcttrunc3", 4.99 );
	set_parameter( "detailed_exceptions", 1 );
	set_parameter( "candidate_blastopts", "-p blastx -e 1E-5 -v 0 -b 1 -F \"\"" );
	set_parameter( "candidate_overlap_ignorestrand",   1 );
	
	return;
}

sub get_mature_peptides {
	my ( $gene ) = @_;

	if ( $$gene{gene_definition} !~ /polyprotein/i ) { return undef }

	my $cleavage_sites;

	$$cleavage_sites{"RR"}{bonus} = 16;
	$$cleavage_sites{"RR"}{offset} = 2;

	$$cleavage_sites{"GG"}{bonus} = 8;
	$$cleavage_sites{"GG"}{offset} = 2;
		
	$$cleavage_sites{"VV"}{bonus} = 8;
	$$cleavage_sites{"VV"}{offset} = 2;
		
	$$cleavage_sites{"AA"}{bonus} = 8;
	$$cleavage_sites{"AA"}{offset} = 2;
		
	$$cleavage_sites{"YS"}{bonus} = 12;
	$$cleavage_sites{"YS"}{offset} = 2;
		
	$$cleavage_sites{"QR"}{bonus} = 8;
	$$cleavage_sites{"QR"}{offset} = 2;

	$$cleavage_sites{"TA"}{bonus} = 4;
	$$cleavage_sites{"TA"}{offset} = 2;

	$$cleavage_sites{"GA"}{bonus} = 4;
	$$cleavage_sites{"GA"}{offset} = 2;

	return mapto_polyprotein( $gene,
		"/usr/local/projects/GCV/Shiliang.dir/YellowFeverV.dir/YFV_matPep.fasta",
		$cleavage_sites );
}

# ---- standardize reference data ----
# return standardized reference name (name is NOT required)
sub get_reference_name {
	my ($reference_id) = @_;

	my $seq = get_reference_seq($reference_id);
	if ( ! defined $seq ) { return undef }
	
	my $id = $$seq{defline};
	if ( $id =~ /polyprotein/i ) { return "POL" }
		
	return $reference_id;
}

sub get_reference_definition {
	my ($reference_id) = @_;

	my $seq = get_reference_seq($reference_id);
	if ( ! defined $seq ) { return undef }

	my $definition = $$seq{defline};
	if ( $definition =~ /polyprotein/i ) { return "polyprotein" }
    
    $definition =~ s/\t/ /g;	
	$definition =~ s/  +/ /g;
	$definition =~ s/^[^ ]*  *//;

	return $definition;
}
}
1;
