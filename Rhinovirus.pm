use strict;
{no warnings 'redefine';
our $myData;
sub customize_defaults {

	# default parameter overrides
	set_parameter( "reference_db", "$myData/Rhinovirus_pep_long.fasta" );
	set_parameter( "min_generef_pctcoverage", 95 );
	set_parameter( "min_generef_pctsimilarity", 0 );
	set_parameter( "min_generef_pctidentity", 0 );
	set_parameter( "max_generef_pcttrunc5", 5 );
	set_parameter( "max_generef_pcttrunc3", 4.99 );
	set_parameter( "detailed_exceptions", 1 );
	set_parameter( "candidate_blastopts", "-p blastx -e 1E-5 -v 0 -b 25 -F \"\"" );
	set_parameter( "candidate_overlap_ignorestrand",   1 );
	
	return;
}

sub get_mature_peptides {
	my ( $gene ) = @_;

	if ( $$gene{gene_definition} !~ /polyprotein/i ) { return undef }

	return mapto_polyprotein( $gene );
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
