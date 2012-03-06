use strict;
our $myData;

{no warnings 'redefine';

my $save_rawgene_reblastopts;

# ---- standardize reference data ----
# return standardized reference name (name is NOT required)
sub get_reference_name {
	my ($reference_id) = @_;
	return $reference_id;	
}

# return standardized reference definition (definition IS required)
sub get_reference_definition {
	my ($reference_id) = @_;

	my $seq = get_reference_seq($reference_id);
	my $definition = $$seq{defline};
	$definition =~ s/\s+/ /g;
	$definition =~ s/^ +//;
	$definition =~ s/^[^ ]* *//;
	$definition =~ s/ +$//;	

	return $definition;
}

sub before_find_rawgenes {
	my ( $genome, $candidate ) = @_;
	
#	$save_rawgene_reblastopts = get_parameter( "rawgene_reblastopts" );

	if ( $$candidate{subject_id} eq "ACN88084.1"
			|| $$candidate{subject_id} eq "ACN88084.1"
			|| $$candidate{subject_id} eq "ACN88085.1"
			|| $$candidate{subject_id} eq "ACN88086.1"
			|| $$candidate{subject_id} eq "ACN88089.1"
			|| $$candidate{subject_id} eq "ACN88091.1"
			|| $$candidate{subject_id} eq "ACN88093.1"
			|| $$candidate{subject_id} eq "ACN88094.1"
			|| $$candidate{subject_id} eq "ACN88100.1"
			|| $$candidate{subject_id} eq "ACN88116.1"
			|| $$candidate{subject_id} eq "ACN88127.1"
			|| $$candidate{subject_id} eq "ACN88132.1"
			|| $$candidate{subject_id} eq "ACN88134.1"
			|| $$candidate{subject_id} eq "ACN88135.1"
			|| $$candidate{subject_id} eq "ACN88136.1"
			|| $$candidate{subject_id} eq "ACN88137.1" ) {
		set_parameter( "allow_splicing", 1 );
	}
	else {
		set_parameter( "allow_splicing", 0 );
	}

#	if (  get_parameter( "allow_splicing" ) ) {	
#		if ( @{$$candidate{hsps}} < 2 ) {
#			set_parameter( "rawgene_reblastopts" , $save_rawgene_reblastopts . " -g F" );
#		}
#	}
}

#sub after_find_rawgenes {
#	my ( $genome, $candidate, $rawgenes ) = @_;
#	
#	set_parameter( "rawgene_reblastopts" , $save_rawgene_reblastopts );
#	
#	if ( ! defined $rawgenes ) { return }
#
#	my @keep;
#	for my $rawgene ( @$rawgenes ) {
##		if ( $$rawgene{subject_id} eq "ACN88088.1" ) { next }
##		if ( $$rawgene{subject_id} eq "ACN88089.1" ) { next }
##		if ( $$rawgene{subject_id} eq "ACN88091.1" ) { next }
##		if ( $$rawgene{subject_id} eq "ACN88095.1" ) { next }
##		if ( $$rawgene{subject_id} eq "ACN88096.1" ) { next }
##		if ( $$rawgene{subject_id} eq "ACN88097.1" ) { next }
##		if ( $$rawgene{subject_id} eq "ACN88099.1" ) { next }
##		if ( $$rawgene{subject_id} eq "ACN88101.1" ) { next }
#		if ( $$rawgene{subject_id} eq "ACN88103.1" && $$rawgene{subject_left} >= 10 ) { next }
##		if ( $$rawgene{subject_id} eq "ACN88132.1" ) { next }
##		if ( $$rawgene{subject_id} eq "ACN88134.1" ) { next }
#		push @keep, $rawgene;
#	}
#	@$rawgenes = @keep;
#}

sub screen_rawgene_variants {
	my ( $genome, $alternates ) = @_;
	
	my $alt;
	for my $alternate ( @$alternates ) {
		if ( $$alternate{subject_id} ne "ACN88094.1" 
				&& $$alternate{subject_id} ne "ACN88100.1" ) { return }
		if ( $$alternate{query_right} > 13800 && $$alternate{query_right} < 14400 ) { $alt = $alternate }			
	}
	if ( defined $alt ) {
		@$alternates = ( $alt );
	}
}
}
1;
