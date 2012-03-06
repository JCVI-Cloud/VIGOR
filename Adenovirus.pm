use strict;
{no warnings 'redefine';
	
my $save_rawgene_reblastopts;

# ---- standardize reference data ----
# return standardized reference name (name is NOT required)
sub screen_rawgene_permutations {
	my ( $genome, $alternates ) = @_;
	
	if ( ! @$alternates ) { return }
	if ( $$alternates[0]{subject_id} ne "ACN88094.1"
			&& $$alternates[0]{subject_id} ne "ACN88100.1"
			&& $$alternates[0]{subject_id} ne "ACN88093.1" ) {
		return;
	}

	my @tmp;
	for my $alternate ( @$alternates ) {
		if ( $$alternate{subject_id} eq "ACN88094.1" || $$alternate{subject_id} eq "ACN88100.1" ) {
			if ( $$alternate{query_right}+$$genome{frag_offset} > 13800 && $$alternate{query_right}+$$genome{frag_offset} < 14800 ) {
				push @tmp, $alternate;
			}
		}
		elsif ( $$alternate{subject_id} eq "ACN88093.1" ) {
			if ( $$alternate{query_right}+$$genome{frag_offset} > 5000 && $$alternate{query_right}+$$genome{frag_offset} < 6000 ) {
				push @tmp, $alternate;
			}
		}	
	}
	if ( @tmp ) {
		@$alternates = @tmp;
	}
}
}
1;
