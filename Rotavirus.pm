use strict;
our $myData;

{no warnings 'redefine';
my $hi_rawgene_sbjcoverage = 80;
my $hi_rawgene_pctsimilarity = 80;

sub customize_defaults {

	# default parameter overrides
	set_parameter( "reference_db", "$myData/rotaVpepDB" );
	set_parameter( "min_generef_pctcoverage", 95 );
	set_parameter( "min_generef_pctsimilarity", 0 );
	set_parameter( "min_generef_pctsimilarity25", 0 );
	set_parameter( "min_generef_pctidentity", 0 );
	set_parameter( "max_generef_pcttrunc5", 5 );
	set_parameter( "max_generef_pcttrunc3", 4.99 );
	set_parameter( "detailed_exceptions", 1 );
	set_parameter( "candidate_blastopts", "-p blastx -e 1E-5 -v 0 -b 25 -F \"\"" );
#	set_parameter( "ignore_strand",   1 );
	
	return;
}

# customized selection of raw genes to be polished into genes
# (default it best bit score)
#sub select_rawgenes_exception {
#	my ( $genome, $candidates, $possible_genes ) = @_;
#
#	my @rawgenes = sort { $$b{bit_score} <=> $$a{bit_score} } @$possible_genes;
#	my $frameSelection;
#	for my $rawgene ( @rawgenes ) {
#		my $frame = $$rawgene{query_frame}; 
#
#		# use this hit?
#		if ( ! defined $frameSelection ) {
#			$$frameSelection{$frame} = $rawgene;
#		}
#		elsif ( exists $$frameSelection{$frame} ) {
#			if ( $$rawgene{vigor_matchwgt} >
#				$$frameSelection{$frame}{vigor_matchwgt} )
#			{
#				$$frameSelection{$frame} = $rawgene;
#			}
#		}
#		elsif ($$rawgene{pct_similarity} >= $hi_rawgene_pctsimilarity
#			&& $$rawgene{pct_scoverage} >= $hi_rawgene_sbjcoverage )
#		{
#			$$frameSelection{$frame} = $rawgene;
#		}
#	}
#
#	my @selection = values %$frameSelection;
#	return \@selection;
#}


# customized post-processing of gene
# add serotype
sub after_find_gene {
	my ( $genome, $rawgene, $gene ) = @_;

	my $genome_seq = $$genome{sequence};
	
	my %serotypes;
	$serotypes{VP7}{db} = "$myData/rotaV_VP7.fasta";
	$serotypes{VP7}{min_identity} = 89;
	$serotypes{VP7}{regexp}       = "\\sVP7(\\.G[0-9]+)\\s";
	$serotypes{VP7}{template}     = "<seroval>";

	$serotypes{VP4}{db} = "$myData/rota_VP4pepDB.fasta";
	$serotypes{VP4}{min_identity} = 89;
	$serotypes{VP4}{regexp}       = "\\sP([0-9]+)\\s";
	$serotypes{VP4}{template}     = ".P[<seroval>]";

	my $blastopts             = "-p blastp -e 1e-10 -F \"\" -v 1 -b 1";
	my $get_alignment_strings = 0;

	for my $serotype ( keys %serotypes ) {
		if ( $$gene{gene_name} =~ /$serotype/i ) {
			my @hits = blast_sequence(
				$$gene{gene_id},           $$gene{protein},
				$serotypes{$serotype}{db}, $blastopts,
				$get_alignment_strings
			);

			#print "\nSEROTYPE $serotype\n";
			#print_hits( @hits );
			for my $hit (@hits) {
				if ( $$hit{pct_identity} > $serotypes{$serotype}{min_identity} )
				{
					if ( $$hit{subject_definition} =~
						/$serotypes{$serotype}{regexp}/i )
					{

						#print "  defline: $$hit{subject_definition}\n"
						#	. "  regexp: $serotypes{$serotype}{regexp}\n"
						#	. "  template: $serotypes{$serotype}{template}\n"
						#	. "  value: $1\n";
						my $value   = $1;
						my $seroval = $serotypes{$serotype}{template};
						$seroval =~ s/<seroval>/$value/;

						#print "  seroval: $seroval\n";
						$$gene{gene_name} .= $seroval;
						last;
					}
				}
			}
		}
	}
	return;
}

# ---- standardize reference data ----
# return standardized reference name (name is NOT required)
sub get_reference_name {
	my ($reference_id) = @_;

	my $acc;

	# default to first word of defline, after the identifier
	my $seq = get_reference_seq($reference_id);
	if ( ! defined $seq ) { return $reference_id }

	my $defline = $$seq{defline};
	$defline =~ s/\t/ /g;
	my @tmp = split /  */, $defline;
	shift @tmp;
	$acc = shift @tmp;

	if ( $defline =~ /^.+\(([NV]\w+)\)/ ) {
		$acc = $1;
	}
	elsif ( $defline =~ /^.+(NSP\d\-\d)/ ) {
		$acc = $1;
	}
	elsif ( $defline =~ /^.+\s+(NS\w+)/ ) {
		$acc = $1;
	}
	elsif ( $defline =~ /^.+\s+(VP\d+)/ ) {
		$acc = $1;
	}

	return $acc;
}

sub get_reference_definition {
	my ($reference_id) = @_;

	my $seq = get_reference_seq($reference_id);
	if ( ! defined $seq ) { return "hypothetical protein" }

	my $definition = $$seq{defline};
	$definition =~ s/\t/ /g;
	$definition =~ s/   */ /g;
	my $id = $definition;
	$id =~ s/ .*$//;

	if ( " $definition " =~ /\.pep/ ) {
		$definition =~ s/^[^ ]*  *//;
		$definition =~ s/^ *$id  *//;
	}
	elsif ( " $definition " =~ /\WVP6\W/i ) {
		$definition = "Inner Capsid Protein VP6, Species\-specific Antigen";
	}
	elsif ( " $definition " =~ /\WVP1\W/i ) {
		$definition = "RNA\-dependent RNA Polymerase";
	}
	elsif ( " $definition " =~ /\WVP2\W/i ) {
		$definition = "RNA Viral Genome Binding Protein VP2";
	}
	elsif ( " $definition " =~ /\WVP3\W/i ) {
		$definition = "Guanylyl Transferase\, mRNA Capping Enzyme";
	}
	elsif ( " $definition " =~ /\WVP4\W/i ) {
		$definition = "Outer Capsid Protein VP4, Virulence Protein";
	}
	elsif ( " $definition " =~ /\WVP7\W/i ) {
		$definition = "Glycoprotein VP7, Neutralization Antigen";
	}
	elsif ( " $definition " =~ /\WNSP1\W/i ) {
		$definition = "Non\-structural Protein NSP1, RNA Binding Protein";
	}
	elsif ( " $definition " =~ /\WNSP2\W/i ) {
		$definition = "Non\-structural Protein NSP2, NTPase";
	}
	elsif ( " $definition " =~ /\WNSP3\W/i ) {
		$definition = "Non\-structural Protein NSP3 binding to mRNA";
	}
	elsif ( " $definition " =~ /\WNSP4\W/i ) {
		$definition = "Non\-structural Protein NSP4, Enterotoxin";
	}
	elsif ( " $definition " =~ /\WNSP5\W/i ) {
		$definition = "Non\-structural Protein NSP5";
	}
	elsif ( " $definition " =~ /\WNSP6\W/i ) {
		$definition =
		  "Non\-structural Protein NSP6, Nucleic Acid Binding Protein";
	}
	else {
		$definition =~ s/^[^ ]*  *//;
		$definition =~ s/^ *$id  *//;
	}

	return $definition;
}
}
1;
