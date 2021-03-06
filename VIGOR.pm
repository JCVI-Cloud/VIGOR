use strict;
our $myBin;
our $myData;
our $myConf;

########################################
# globals
#my $rootspace = "/usr/local/scratch/VIRAL/VIGOR/tempspace";
my $rootspace = "/usr/local/scratch/vigor/tempspace";

#my $rootspace = "/home/jhoover/VIRAL/VIGOR";

my %codons;
my %cleavage_sites;

my %splice_pairs;
my %reference_seqs;
my %refmat;
my %parameters;
my %slippageGenes;
my %splicedGenes;
my %readthruGenes;
my ( $next_sequence_buffer, $next_sequence_eof );
my ( @sequence_fragments, $sequence_fragcount, $frag_sequence_id );
my $profileid = 0;

sub initialize_defaults {
	########################################
	# use default codons/splice sites
	set_default_codons();
	set_default_splice_pairs();
	set_default_cleavage_sites();

	########################################
	# set parameter defaults
	my $vigorspace = create_workspace();
	set_parameter( "vigorspace", $vigorspace);

	# general
	set_parameter( "allow_partial_genes", 1 );
	set_parameter( "expected_genome_length", 12000000000 );
	set_parameter( "ignore_strand", 0 );		# 0=consider strand, 1=ignore
	set_parameter( "max_aa_gap", 10 );

	# frameshifts
	set_parameter( "frameshift_detection", 1 );
	set_parameter( "allow_ribosomal_slippage", 2 );	# 0=do not allow, 1=allow, 2=check deflin/exceptions
	set_parameter( "ribosomal_slippage_list", "" );

	# splicing
	set_parameter( "allow_splicing", 2 );		# 0=do not allow, 1=allow, 2=check deflin/exceptions
	set_parameter( "spliced_gene_list", "" );	# comma separated list of spliced genes (overrides allow_splicing),
												# values should match values returned by get_reference_name()
	set_parameter( "min_intron_size", 100 );	# relevant only when splicing allowed
	set_parameter( "max_intron_size", 7000 );	# relevant only when splicing allowed
	set_parameter( "min_exon_size", 3 );		# relevant only when splicing allowed
	set_parameter( "small_exon_autofind", 1 );
	set_parameter( "splicesite_search_range", 100 );

	# candidate regions
	set_parameter( "tile_size", 100000 );
	set_parameter( "tile_overlap", 10000 );
	set_parameter( "hspmerge_blastopts", "-p blastx -e 0.01 -W 3 -P 1 -F \"\" -b 1 -v 0 " );

	set_parameter( "reference_db", "$myData/viral_refpep" );
	set_parameter( "candidate_refdb", "" );
	set_parameter( "cluster_candidaterefs", 0 );		# cluster candidate ref sequences
	set_parameter( "candidateref_clusterL", 0.90 );		# mutual coverage
	set_parameter( "candidateref_clusterS", 80 );		# %identity
	set_parameter( "candidate_blastopts", "-p blastx -e 1E-10 -v 0 -F \"\"" );
	set_parameter( "min_candidate_pctsimilarity", 40 );
	set_parameter( "min_candidate_sbjcoverage", 50 );
	set_parameter( "candidate_overlap_method", 1 );			# 0=span, 1=hsp (equivalent unless splicing allowed)
	set_parameter( "candidate_overlap_threshold", 0.80 );	# 80% of query range
	set_parameter( "candidate_overlap_ismutual", 0 );		# 0=either
															# 1=both seqs meet threshold unless same gene
															# 2=both seqs meet threshold

	# rawgenes
	set_parameter( "rawgene_reblastopts", "-p blastx -e 1e-10 -b 20 -v 0 -F \"\"" );
	set_parameter( "rawgene_extend_span", 50 );
	set_parameter( "rawgene_extend_hsp", 25 );
	set_parameter( "min_rawgene_pctsimilarity", 50 );
	set_parameter( "min_rawgene_sbjcoverage", 50 );
	set_parameter( "rawgene_overlap_threshold", 0.80 );	# 80% of query range
	set_parameter( "rawgene_overlap_ismutual", 0 );		# 0=either, 1=both seqs meet threshold

	# genes
	set_parameter( "allow_fuzzy_matpeps", 1 );			# allow "fuzzy" (<>) edges on matpeps, 1=yes, 0=no
	set_parameter( "startcodon_search_aarange", "50%:25%" );
	set_parameter( "stopcodon_search_aaextension", "100%" );
	set_parameter( "allow_stopcodon_readthru", 0 );	# 0=do not allow, 1=allow, 2=check deflin/exceptions
	set_parameter( "stop_readthru_list", "" );
	set_parameter( "pep_bl2seqopts", "-p blastx -F \"\" -e 1e-5 -d 3000000000" );

	# polyproteins
	set_parameter( "polyprotein_threshold", 2500 );		# proteins >= are automatically annotated for matpeps
	set_parameter( "mature_pep_refdb", "$myData/viruses.mature_peps" );
	set_parameter( "mature_pep_blastopts", "-p blastp -e 100 -b 4000 -v 0" ); # -F \"\"" );
	set_parameter( "mature_pep_mincoverage", 25 );
	set_parameter( "mature_pep_minsimilarity", 25 );

	# validation
	set_parameter( "min_generef_pctcoverage", 66 );
	set_parameter( "min_generef_pctsimilarity", 55 );
	set_parameter( "min_generef_pctidentity", 40 );
	set_parameter( "max_generef_pcttrunc5", 5 );
	set_parameter( "max_generef_pcttrunc3", 33 );
	set_parameter( "min_generef_pctsimilarity25", 55 );
	set_parameter( "detailed_exceptions", 1 );
}

########################################
# manipulate parameter set
sub read_config_file {
	my ( $cfgfile ) = @_;

	open( CFG, "<$cfgfile" ) || die "Could not read configuration file \"cfgfile\"\n";

	while ( my $line = <CFG> ) {
		chomp $line;
		if ( $line =~ /^ *#/ ) { next }
		$line =~ s/^ *//;
		$line =~ s/ *$//;
		if ( ! $line ) { next }
		my ( $name, @values ) = split /=/, $line;
		if ( defined $name ) { set_parameter( $name, join( "=", @values ) ) }
	}
	close CFG;
}

sub set_parameter {
	my ( $name, $value ) = @_;
	$parameters{lc $name} = $value;
	return;
}

sub get_parameter {
	my ($name) = @_;
	return $parameters{lc $name};
}

sub show_parameters {
	print "\nparameters\n";
	for my $param ( sort keys %parameters ) {
		print "\t$param\t$parameters{$param}\n";
	}
	print "\n";
	return;
}

sub get_parameters {
	return \%parameters;
}

sub set_parameters {
	my ($new_parameters) = @_;
	%parameters = %$new_parameters;
	return;
}

########################################
# manipulate codon set
sub set_default_codons {

	$codons{'GC.'}         = 'A';    # Alanine
	$codons{'TG[TCY]'}     = 'C';    # Cysteine
	$codons{'GA[TCY]'}     = 'D';    # Aspartic Acid
	$codons{'GA[AGR]'}     = 'E';    # Glutamic Acid
	$codons{'TT[TCY]'}     = 'F';    # Phenylalanine
	$codons{'GG.'}         = 'G';    # Glycine
	$codons{'CA[TCY]'}     = 'H';    # Histidine
	$codons{'AT[TCAY]'}    = 'I';    # Isoleucine
	$codons{'AA[AGR]'}     = 'K';    # Lysine
	$codons{'TT[AGR]|CT.'} = 'L';    # Leucine
	$codons{'ATG'}         = 'M';    # Methionine
	$codons{'AA[TCY]'}     = 'N';    # Asparagine
	$codons{'CC.'}         = 'P';    # Proline
	$codons{'CA[AGR]'}     = 'Q';    # Glutamine
	$codons{'CG.|AG[AGR]'} = 'R';    # Arginine
	$codons{'TC.|AG[TCY]'} = 'S';    # Serine
	$codons{'AC.'}         = 'T';    # Threonine
	$codons{'GT.'}         = 'V';    # Valine
	$codons{'TGG'}         = 'W';    # Tryptophan
	$codons{'TA[TCY]'}     = 'Y';    # Tyrosine
	$codons{'TA[AGR]|TGA'} = '*';    # Stop
}

sub get_aa_codons {

	my %tmp;
	for my $regexp ( keys %codons ) {
		my $aa = $codons{$regexp};
		if ( exists $tmp{$aa} ) {
			$tmp{$aa} .= "|$regexp";
		}
		else {
			$tmp{$aa} = $regexp;
		}
	}

	my %aacodons;
	for my $aa ( keys %tmp ) {
		my $regexp = $tmp{$aa};
		my $dna = regexp2dna( $regexp );
		$aacodons{$aa} = $dna;
		print "aa: $aa\tregexp: $regexp\tdna: $dna\n";
	}

	return %aacodons;
}

sub regexp2dna {
	my ( $regexp ) = @_;

	my $dna = $regexp;
	$dna =~ s/\./N/g;
	$dna =~ s/\[[^]]+\]/N/g;
	my @tmp = split /\|/, $dna;
	$dna = "";
	for my $i ( 0..2 ) {
		my $na;
		for my $tmp( @tmp ) {
			if ( !defined $na ) {
				$na = substr( $tmp, $i, 1 );
			}
			elsif ( $na ne substr( $tmp, $i, 1 ) ) {
				$na = "N";
				last;
			}
		}
		$dna .= $na;
	}
	return $dna;
}

sub add_codon {
	my ( $codon_regexp, $aa ) = @_;

	$codons{$$codon_regexp} = $aa;
	return;
}

sub get_codons {
	return \%codons;
}

sub set_codons {
	my ($new_codons) = @_;
	%codons = %$new_codons;
	return;
}

sub set_default_cleavage_sites {
	$cleavage_sites{"GG"}{bonus} = 10;
	$cleavage_sites{"GG"}{offset} = 2;
	$cleavage_sites{"RR"}{bonus} = 7;
	$cleavage_sites{"RR"}{offset} = 2;
	$cleavage_sites{"QA"}{bonus} = 10;
	$cleavage_sites{"QA"}{offset} = 1;
	$cleavage_sites{"QS"}{bonus} = 10;
	$cleavage_sites{"QS"}{offset} = 1;
	$cleavage_sites{"QN"}{bonus} = 5;
	$cleavage_sites{"QN"}{offset} = 1;
}

sub get_cleavage_sites {
	return \%cleavage_sites;
}


sub set_cleavage_sites {
	my ($new_sites) = @_;
	%cleavage_sites = %$new_sites;
	return;
}

########################################
# manipulate splice donor/acceptor sets
sub set_default_splice_pairs {
	$splice_pairs{GT}{AG} = 0;
	$splice_pairs{GC}{AG} = 0.25;
	$splice_pairs{AT}{AC} = 0.50;
	$splice_pairs{GT}{CA} = 0.75;
	return;
}

sub get_splice_pairs {
	return \%splice_pairs;
}

sub set_splice_pairs {
	my ($new_splice_pairs) = @_;
	%splice_pairs = %$new_splice_pairs;
	return;
}

sub add_splice_pair {
	my ( $penalty, $donor, @acceptors ) = @_;

	my $d = uc $donor;
	for my $a ( @acceptors ) {
		$splice_pairs{$d}{uc $a} = $penalty;
	}
}

sub get_splice_donors {
	return \%splice_pairs;
}

sub get_splice_acceptors {
	my %acceptors;
	for my $d ( keys %splice_pairs ) {
		my %acctmp = %{$splice_pairs{$d}};
		for my $a ( keys %acctmp ) {
			$acceptors{$a}{$d} = $splice_pairs{$d}{$a};
		}
	}
	return \%acceptors;
}

########################################
# reference sequences
sub get_reference_seq {
	my ($seq_id) = @_;
	return $reference_seqs{$seq_id};
}

sub add_reference_seq {
	my ( $ref ) = @_;
	$reference_seqs{$$ref{id}} = $ref;
}

sub get_reference_seqs {
	return \%reference_seqs;
}

sub get_reference_dbsize {
	return get_parameter( "reference_dbsize" );
}

sub get_reference_db {
	return get_parameter( "reference_db" );
}

sub load_reference_db {

	my $reference_db = get_parameter( "reference_db" );
	my $vigorspace = get_parameter("vigorspace");

	# convert gbk/fasta to blast database
	my $ref_fasta  = "$vigorspace/ref.fasta";
	if ( $reference_db =~ /\.gbk/ ) {
#print "$myBin/gbk_to_faa.pl < $reference_db > $ref_fasta";
		system "cat $reference_db | $myBin/gbk_to_faa.pl > $ref_fasta";
#print "formatdb -i $ref_fasta -p T";
		system "formatdb -i $ref_fasta -p T -l $vigorspace/formatdb.log";
		set_parameter( "reference_db", $ref_fasta );
		$reference_db = $ref_fasta;
#print `grep ">" $ref_fasta` . "\n";
	}
	elsif ( $reference_db =~ /\.fa(st){0,1}a/ ) {
		system "cp $reference_db $ref_fasta";
		system "formatdb -i $ref_fasta -p T -l $vigorspace/formatdb.log";
		set_parameter( "reference_db", $ref_fasta );
		$reference_db = $ref_fasta;
	}

	# load reference sequences
	else {
		system "fastacmd -d $reference_db -D 1 | sed 's/\t/ /g' | sed 's/^>[^ ]* */>/' > $ref_fasta";
	}

	%reference_seqs   = loadFasta($ref_fasta);
	my $reference_dbsize = get_db_size($reference_db);
	set_parameter( "reference_dbsize", $reference_dbsize );

	# has a special reference db been specified for candidate region identification?
	my $candidate_db = get_parameter( "candidate_refdb" );

	# cluster standard reference set for candidate region identification
	my $cluster = get_parameter( "cluster_candidaterefs" );
#print "refdb=$reference_db  candb=$candidate_db\n";
	if ( $cluster ) {
		my $L = get_parameter( "candidateref_clusterL" );
		my $S = get_parameter( "candidateref_clusterS" );
		$candidate_db = "$vigorspace/candidates";
		fasta_to_clusterdb( $ref_fasta, $candidate_db, $L, $S );
	}

	# use default reference db for candidate and rawgene identification
	elsif ( $candidate_db le " " ) {
		$candidate_db = $reference_db;
	}
	set_parameter( "candidate_refdb", $candidate_db );
	my $candidate_dbsize = get_db_size($candidate_db);
	set_parameter( "candidate_dbsize", $candidate_dbsize );
#print "refdb=$reference_db  candb=$candidate_db\n";
}

sub fasta_to_clusterdb {
	my ( $infasta, $cldbpath, $cluster_L, $cluster_S ) = @_;

	my $vigorspace = get_parameter("vigorspace");

	my $cmd = "blastclust -i $infasta -S $cluster_S -L $cluster_L | sed 's/ .*\$//'";
#print "$cmd\n";
	my @clusters = split /\n/, `$cmd`;
	my %refs = loadFasta( $infasta );
	open( DB, ">$cldbpath" );
	for my $cluster ( @clusters ) {
#print "$cluster\n";
		if ( exists $refs{$cluster} ) {
			print DB ">$refs{$cluster}{defline}\n$refs{$cluster}{sequence}\n";
		}
	}
	close DB;
	$cmd = "formatdb -i $cldbpath -p T -l $cldbpath.log";
	system $cmd;
}

sub blastdb_to_clusterdb {
	my ( $indb, $cldbpath, $cluster_L, $cluster_S ) = @_;

	my $vigorspace = get_parameter("vigorspace");

	my $infasta  = "$vigorspace/tmpdb.fasta";
	system "fastacmd -d $indb -D 1 | sed 's/\t/ /g' | sed 's/^>[^ ]* */>/' > $infasta";

	fasta_to_clusterdb( $infasta, $cldbpath, $cluster_L, $cluster_S );
	unlink $infasta;
}

########################################
# miscellaneous utilities
sub minval {
	my (@vals) = @_;

	if ( !@vals ) { return undef }
	my $val = $vals[0];
	for my $i ( 1 .. @vals - 1 ) {
		if ( $vals[$i] < $val ) { $val = $vals[$i] }
	}

	return $val;
}

sub maxval {
	my (@vals) = @_;

	if ( !@vals ) { return undef }
	my $val = $vals[0];
	for my $i ( 1 .. @vals - 1 ) {
		if ( $vals[$i] > $val ) { $val = $vals[$i] }
	}

	return $val;
}

sub create_workspace {
	my $vigorspace;
	if ( !-e $rootspace ) { mkdir $rootspace }

	$vigorspace = $rootspace . "/" . rand(1000000);

#print "vigorspace=$vigorspace\n";
	while ( -e $vigorspace ) {
		$vigorspace = $rootspace . "/" . rand(1000000);
#print "vigorspace=$vigorspace\n";
	}

	mkdir $vigorspace;
	return $vigorspace;
}

########################################
# various comparators used in sorts
sub compare_hsp_seqs {
	my ( $a, $b ) = @_;

	if ( $$a{query_id} lt $$b{query_id} ) {
		return -1;
	}
	elsif ( $$a{query_id} gt $$b{query_id} ) {
		return 1;
	}
	elsif ( $$a{subject_id} lt $$b{subject_id} ) {
		return -1;
	}
	elsif ( $$a{subject_id} gt $$b{subject_id} ) {
		return 1;
	}
	elsif ( $$a{orientation} < $$b{orientation} ) {
		return -1;
	}
	elsif ( $$a{orientation} > $$b{orientation} ) {
		return 1;
	}
	else {
		return 0;
	}
}

sub conflict_compare {
	my ( $a, $b ) = @_;

	if ( $$a{query_id} lt $$b{query_id} ) {
		return -1;
	}
	elsif ( $$a{query_id} gt $$b{query_id} ) {
		return 1;
	}
	elsif ( $$a{orientation} > $$b{orientation} ) {
		return -1;
	}
	elsif ( $$a{orientation} < $$b{orientation} ) {
		return 1;
	}
	elsif ( $$a{subject_id} lt $$b{subject_id} ) {
		return -1;
	}
	elsif ( $$a{subject_id} gt $$b{subject_id} ) {
		return 1;
	}
	elsif ( $$a{bit_score} > $$b{bit_score} ) {
		return -1;
	}
	elsif ( $$a{bit_score} < $$b{bit_score} ) {
		return 1;
	}
	elsif ( $$a{position_penalty} < $$b{position_penalty} ) {
		return -1;
	}
	elsif ( $$a{position_penalty} > $$b{position_penalty} ) {
		return 1;
	}
	else {
		return 0;
	}
}

sub order_by_seqs_and_scores {
	my ( $a, $b ) = @_;

	my $seqcmp = compare_hsp_seqs( $a, $b );
	if ( $seqcmp != 0 ) { return $seqcmp }

	if ( $$a{bit_score} > $$b{bit_score} ) {
		return -1;
	}
	elsif ( $$a{bit_score} < $$b{bit_score} ) {
		return 1;
	}
	elsif ( $$a{evalue} < $$b{evalue} ) {
		return -1;
	}
	elsif ( $$a{evalue} > $$b{evalue} ) {
		return 1;
	}
	else {
		return 0;
	}
}

sub order_by_seqs_and_positions {
	my ( $a, $b ) = @_;

	my $cmp = compare_hsp_seqs( $a, $b );
	if ( $cmp != 0 ) { return $cmp }

	$cmp = compare_qry_positions( $a, $b );
	if ( $cmp != 0 ) { return $cmp }

	return compare_sbj_positions( $a, $b );
}

sub compare_qry_positions {
	my ( $a, $b ) = @_;

	if ( $$a{orientation} > $$b{orientation} ) {
		return 1;
	}
	elsif ( $$a{orientation} < $$b{orientation} ) {
		return -1;
	}

	my $astart = $$a{query_left};
	my $aend = $$a{query_right};
	my $bstart = $$b{query_left};
	my $bend = $$b{query_right};

	if ( $$a{orientation} == -1 ) {
		$astart = -$$a{query_right};
		$aend = -$$a{query_left};
		$bstart = -$$b{query_right};
		$bend = -$$b{query_left};
	}

	if ( $astart  < $bstart ) {
		return -1;
	}
	elsif ( $astart  > $bstart ) {
		return 1;
	}
	elsif ( $aend < $bend ) {
		return -1;
	}
	elsif ( $aend > $bend ) {
		return -1;
	}
	else {
		return 0;
	}
}

sub order_left_to_right {
	my ( $a, $b ) = @_;
	if ( $$a{query_id} lt $$b{query_id} ) {
		return -1;
	}
	elsif ( $$a{query_id} gt $$b{query_id} ) {
		return 1;
	}
	elsif ( $$a{query_left} < $$b{query_left} )
	{
		return -1;
	}
	elsif ( $$a{query_left} > $$b{query_left} )
	{
		return 1;
	}
	elsif ( $$a{orientation} > $$b{orientation} ) {
		return -1;
	}
	elsif ( $$a{orientation} < $$b{orientation} ) {
		return 1;
	}
	elsif ( $$a{query_right} < $$b{query_right} )
	{
		return -1;
	}
	elsif ( $$a{query_right} > $$b{query_right} )
	{
		return 1;
	}
	else {
		return 0;
	}
}

sub compare_gene_ids {
	my ( $a, $b ) = @_;

	my ( @atmp ) = split /\./, $$a{gene_id};
	my $agene_num = pop @atmp;
	my $agenome_id = join( ".", @atmp );
	my ( @btmp ) = split /\./, $$b{gene_id};
	my $bgene_num = pop @btmp;
	my $bgenome_id = join( ".", @btmp );

	my $asuffix = " ";
	my $bsuffix = " ";
	if ( $agene_num =~ /[a-z]/i ) {
		$asuffix = $agene_num;
		$asuffix =~ s/[^a-z]//gi;
		$agene_num =~ s/[a-z]//gi;
	}
	if ( $bgene_num =~ /[a-z]/i ) {
		$bsuffix = $bgene_num;
		$bsuffix =~ s/[^a-z]//gi;
		$bgene_num =~ s/[a-z]//gi;
	}

	if ( $agenome_id ne $bgenome_id ) { return $agenome_id cmp $bgenome_id }
	elsif ( $agene_num != $bgene_num ) { return $agene_num <=> $bgene_num }
	else { return $asuffix cmp $bsuffix }
}

sub compare_gene_positions {
	my ( $a, $b ) = @_;

	if ( $$a{orientation} < $$b{orientation} ) {
		return -1;
	}
	elsif ( $$a{orientation} > $$b{orientation} ) {
		return 1;
	}
	elsif ( $$a{orientation} * $$a{start_codon} <
		$$b{orientation} * $$b{start_codon} )
	{
		return -1;
	}
	elsif ( $$a{orientation} * $$a{start_codon} >
		$$b{orientation} * $$b{start_codon} )
	{
		return 1;
	}
	elsif (
		$$a{orientation} * $$a{stop_site} < $$b{orientation} * $$b{stop_site} )
	{
		return -1;
	}
	elsif (
		$$a{orientation} * $$a{stop_site} > $$b{orientation} * $$b{stop_site} )
	{
		return 1;
	}
	else {
		return 0;
	}
}

sub compare_sbj_positions {
	my ( $a, $b ) = @_;
	if ( $$a{subject_left} < $$b{subject_left} ) {
		return -1;
	}
	elsif ( $$a{subject_left} > $$b{subject_left} ) {
		return 1;
	}
	elsif ( $$a{subject_right} < $$b{subject_right} ) {
		return -1;
	}
	elsif ( $$a{subject_right} > $$b{subject_right} ) {
		return 1;
	}
	else {
		return 0;
	}
}

sub compare_hsps {
	my ( $a, $b ) = @_;
	if ( $$a{subject_left} < $$b{subject_left} ) {
		return -1;
	}
	elsif ( $$a{subject_left} > $$b{subject_left} ) {
		return 1;
	}
	elsif ( $$a{subject_right} > $$b{subject_right} ) {
		return -1;
	}
	elsif ( $$a{subject_right} < $$b{subject_right} ) {
		return 1;
	}
	else {
		return 0;
	}
}

########################################
# sequences and fasta files

# extract subsequence
# reverse complements if start>end
sub subsequence {
	my ( $sequence, $start, $end ) = @_;

	if ( $start < 1 ) { return "" }
	if ( $start > length( $sequence ) ) { return "" }

	if ( $start <= $end ) {
		return substr( $sequence, $start - 1, $end - $start + 1 );
	}
	else {
		my $subseq = substr( $sequence, $end - 1, $start - $end + 1 );
		return reverse_complement( $subseq );
	}
}

sub reverse_complement {
	my ( $seq  ) = @_;

	if ( ! defined $seq ) { return "" }

	my $revcomp = reverse $seq;
	$revcomp =~ tr/AGCTMRWSYKVHDBNagctmrwsykvhdbn/TCGAKYWSRMBDHVNtcgakywsrmbdhvn/;

	return $revcomp;
}

# translate DNA sequence to amino acid sequence
sub DNA2AA {

	my ( $seq, $translation_exception )      = @_;

	if ( ! defined $seq ) { return "" }

	my $len        = length($seq);
	if ( ! $len ) { return "" }

	my $start_site = 0;
	my $protein;
	while ( $start_site <= $len - 3 ) {
		my $codon_dna = substr( $seq, $start_site, 3 );
		my $AA = &codon2aa($codon_dna);

		$protein    = "$protein$AA";
		$start_site = $start_site + 3;
	}

	if ( defined $translation_exception ) {
		my $pos = $$translation_exception{position};
		my $aa = $$translation_exception{aa};
		$protein = substr( $protein, 0, $pos-1 ) . $aa . substr( $protein, $pos );
	}
	return $protein;
}

sub AA2DNA {
	my  ( $pep ) = @_;

	my %aacodons = get_aa_codons();
	my $dna = "";
	for my $i ( 0 .. length( $pep ) - 1 ) {
		my $aa = uc substr( $pep, $i, 1 );
		$dna .= $aacodons{$aa};
	}
	return $dna;
}

sub codon2aa {
	my ( $codon ) = @_;

	for my $key ( keys %codons ) {
		if ( $codon =~ /$key/i ) {
			return $codons{$key};
		}
	}

	# probably contain an IUPAC ambiguity code;
	return "X";
}

# read fasta file into hash;
sub loadFasta {
	my ($fasta) = @_;

	my %contents;
	open( FASTA, "<$fasta" ) || die "\nCould not read fasta file \"$fasta\"\n";

	my $sequence = "";
	my $defline;
	while ( my $line = <FASTA> ) {
		chomp $line;
		$line =~ s/\r//g;
		$line =~ s/\t/ /g;
		if ( $line =~ /^ *$/ ) { next }
		if ( $line =~ /^>/ ) {
			if ( length($sequence) ) {
				my %entry;
				$entry{defline}  = $defline;
				$entry{sequence} = $sequence;
				$entry{length} = length( $sequence );
				if ( $sequence =~ /\*$/ ) { $entry{length}-- }
				my ($id) = get_defline_id( $defline );
				$entry{id} = $id;
				$contents{$id} = \%entry;
			}
			$defline = $line;
			$defline =~ s/^ *> *//;
			$sequence = "";
		}
		else {
			$sequence .= $line;
		}
	}
	close FASTA;
	if ( length($sequence) ) {
		my %entry;
		$entry{defline}  = $defline;
		$entry{sequence} = $sequence;
		$entry{length} = length( $sequence );
		if ( $sequence =~ /\*$/ ) { $entry{length}-- }
		my ($id) = get_defline_id( $defline );
		$entry{id} = $id;
		$contents{$id} = \%entry;
	}


	return %contents;
}

# fetch next sequence from fasta file
sub next_sequence {
	if ($next_sequence_eof) { return undef }
	my ($FASTA) = @_;

	my %entry;
	$entry{id}       = "";
	$entry{defline}  = "";
	$entry{sequence} = "";
	$entry{frag_offset} = 0;
	my $full_length = get_parameter( "expected_genome_length" );

	if ( defined $next_sequence_buffer && length($next_sequence_buffer) ) {
		$next_sequence_buffer =~ s/^ *> *//;
		$entry{defline} = $next_sequence_buffer;
		( $entry{id} ) = split / /, $entry{defline};
	}

	while ( $next_sequence_buffer = <$FASTA> ) {
		chomp $next_sequence_buffer;
		$next_sequence_buffer =~ s/\r//g;
		$next_sequence_buffer =~ s/\t/ /g;
		if ( $next_sequence_buffer =~ /^ *$/ ) { next }
		if ( $next_sequence_buffer =~ /^>/ ) {
			my $seqlen = length( $entry{sequence} );
			if ( $seqlen ) {
				$entry{seqlen} = $seqlen;
				my $full_left = $seqlen - $full_length;
				if ( $full_left > 1 ) { $full_left = 1 }
				$entry{full_left} = $full_left;
				my $full_right = $full_length - $seqlen;
				if ( $full_right > $seqlen ) { $full_right = $seqlen }
				$entry{full_right} = $full_right;
				return \%entry;
			}
			$next_sequence_buffer =~ s/^ *> *//;
			$entry{defline}  = $next_sequence_buffer;
			( $entry{id} ) = split / /, $entry{defline};
			$entry{sequence} = "";
		}
		else {
			$entry{sequence} .= $next_sequence_buffer;
		}
	}

	$next_sequence_eof = 1;
	my $seqlen = length( $entry{sequence} );
	if ( $seqlen ) {
		$entry{seqlen} = $seqlen;
		my $full_left = $seqlen - $full_length;
		if ( $full_left > 1 ) { $full_left = 1 }
		$entry{full_left} = $full_left;
		my $full_right = $full_length - $seqlen;
		if ( $full_right > $seqlen ) { $full_right = $seqlen }
		$entry{full_right} = $full_right;
		return \%entry;
	}
	return undef;
}

sub next_frag {
	my ( $sequence ) = @_;

	if ( ! defined $frag_sequence_id || $$sequence{id} ne $frag_sequence_id ) {
		$sequence_fragcount = undef;
		$frag_sequence_id = $$sequence{id};
	}
	$$sequence{frag_offset} = 0;

	if ( ! defined $sequence_fragcount ) {
		my @fragseqs = split /(N{20,})/, $$sequence{sequence};
		$sequence_fragcount = 0;
		my $offset = 0;
		for my $fragseq ( @fragseqs ) {
			if ( $fragseq =~ /[^N]/ ) {
				my %frag = %$sequence;
				$frag{sequence} = $fragseq;
				$frag{seqlen} = length( $frag{sequence} );
				$frag{frag_offset} = $offset;
				$frag{frag_left} = $offset + 1;
				$frag{frag_right} = $offset + $frag{seqlen};
				$frag{full_left} = $$sequence{full_left} - $offset;
				$frag{full_right} = $$sequence{full_right} + $$sequence{seqlen} - $frag{frag_right};
				push @sequence_fragments, \%frag;
				$sequence_fragcount++;
			}
			else {
				my $gapstart = $offset+1;
				my $gapend = $offset + length( $fragseq );
			}
			$offset += length( $fragseq );
		}
	}
	if ( $sequence_fragcount <= 0 ) {
		return undef;
	}
	else {
		$sequence_fragcount--;
		return shift @sequence_fragments;
	}
}

# parse id from defline
sub get_defline_id {
	my ($defline) = @_;

	my ($rawid) = split / /, $defline;
	$rawid =~ s/[,\| ]*$//;

	my @id = split( /\|/, $rawid );
	if ( scalar @id < 2 ) { return $rawid }

	if ( $id[0] eq "gnl" && @id > 2 ) {
		shift @id;
		shift @id;
		$rawid = join( "|", @id );
		if ( @id <= 2 ) { return $rawid }
	}

	my $i = 0;
	while ( $i < @id - 1 ) {
		if (
			index( ".gi.gb.rf.emb.dbj.pir.prf.sp.ref.",
				"." . lc( $id[$i] ) . "." ) >= 0
		  )
		{
			return lc( $id[$i] ) . "|" . $id[ $i + 1 ];
		}
		$i++;
	}
	return $rawid;
}

# return blast db size
sub get_db_size {
	my ($dbpath) = @_;

	my @tmp = split /\n/, `fastacmd -d $dbpath -I`;
	for my $dbsize (@tmp) {
		if ( $dbsize =~ /total letters/ ) {
			$dbsize =~ s/total letters.*//;
			$dbsize =~ s/^.*; //;
			$dbsize =~ s/[^0-9]//g;
			return $dbsize;
		}
	}

	return undef;
}

########################################
# blast
#
# break sequence into overlapping tiles, blast, and merge results back together
sub blast_sequence {
	my ( $defline, $sequence, $database, $options, $get_alignment_strings,
		$tilesize, $tileoverlap, $dbg )
	  = @_;

	my $vigorspace = get_parameter( "vigorspace" );
	if ( ! defined $dbg ) { $dbg = 0 }

	my @hsps;
	my $tmpfile = "$vigorspace/tile.fasta";
	my $xmlfile = "$vigorspace/tile.xml";
	unlink $tmpfile;
	unlink $xmlfile;

	my %uniqhits;
	for my $tile ( split_sequence( $sequence, $tilesize, $tileoverlap ) ) {
		open( TMP, ">$tmpfile" );
		print TMP ">$defline\n$$tile{sequence}";
		close TMP;
		my $blastcmd =
		  "blastall -i $tmpfile -d $database $options -m 7 > $xmlfile";
		system $blastcmd;
#print "$blastcmd\n";
#print "\n" . `cat $tmpfile` . "\n";
		my @tilehsps = offset_adjustment(
			$$tile{offset} + 1, length( $$tile{sequence} ), length($sequence), 1, undef, undef,
			parse_blastxml( $xmlfile, $get_alignment_strings, $dbg ) );

		for my $hsp ( sort { order_by_seqs_and_positions( $a, $b ) } @tilehsps ) {
			my $key = "$$hsp{subject_id}:$$hsp{subject_left}-$$hsp{subject_right}";
			if ( ! exists $uniqhits{$key}
					|| $uniqhits{$key}{vigor_pctsimilarity} <  $$hsp{vigor_pctsimilarity} ) {
				$uniqhits{$key} = $hsp;
			}
		}
	}
	unlink $tmpfile;
	unlink $xmlfile;

	@hsps = sort { order_left_to_right( $a, $b ) } values %uniqhits;

	return @hsps;
}

# break large input sequence into overlapping tiles
sub split_sequence {
	my ( $sequence, $size, $overlap ) = @_;

	if ( !defined $size ) {
		my %tile;
		$tile{offset}   = 0;
		$tile{sequence} = $sequence;
		my @tiles = ( \%tile );
		return @tiles;
	}
	elsif ( $overlap >= $size ) {
		die "\nsplit_sequence: overlap $overlap >= size $size\n";
	}

	my $seqlen = length($sequence);
	if ( $size >= $seqlen ) {
		my %tile;
		$tile{offset}   = 0;
		$tile{sequence} = $sequence;
		my @tiles = ( \%tile );
		return @tiles;
	}

	my %tiles;
	my $offset = 0;
	while ( $offset < $seqlen ) {
		if ( $offset + $size > $seqlen ) {
			$offset = $seqlen - $size;
			if ( exists $tiles{$offset} ) { last }
		}
		my %tile;
		$tile{offset}   = $offset;
		$tile{sequence} = substr( $sequence, $offset, $size );
		$tiles{$offset} = \%tile;
		$offset += ( $size - $overlap );
	}

	return sort { $$a{offset} <=> $$b{offset} } values %tiles;
}

# parse blast result xml
sub parse_blastxml {
	my ( $xmlfile, $get_alignment_strings, $dbg ) = @_;

	if ( ! defined $dbg ) { $dbg = 0 }

	#	my $euler = 2.718281828;
	my %xmlattr;
	$xmlattr{"Iteration_query-def"} = "query_definition";
	$xmlattr{"Iteration_query-len"} = "query_length";
	$xmlattr{"Hit_def"}             = "subject_definition";
	$xmlattr{"Hit_len"}             = "subject_length";
	$xmlattr{"Hsp_bit-score"}       = "bit_score";

	#	$xmlattr{"Hsp_score"} = "hsp_score";
	$xmlattr{"Hsp_evalue"}     = "evalue";
	$xmlattr{"Hsp_query-from"} = "query_left";
	$xmlattr{"Hsp_query-to"}   = "query_right";
	$xmlattr{"Hsp_hit-from"}   = "subject_left";
	$xmlattr{"Hsp_hit-to"}     = "subject_right";

	$xmlattr{"Hsp_query-frame"} = "query_frame";
	$xmlattr{"Hsp_hit-frame"}   = "subject_frame";
	$xmlattr{"Hsp_identity"}    = "num_identical";
	$xmlattr{"Hsp_positive"}    = "num_similar";
	$xmlattr{"Hsp_gaps"}        = "num_gaps";
	$xmlattr{"Hsp_align-len"}   = "alignment_length";
	$xmlattr{"Hsp_qseq"}        = "query_alignstr";
	my $last_attribute = "query_alignstr";

	# bit flags signal which alignment strings to retrieve:
	# 0=none, 1 = query, 2 = subject, 4 = midline
	# (i.e. 3 = query & subject alignment strings)
	#
	# note: we always EXTRACT the query string so we can
	# count stop codons, but we only keep it if flagged
	if ($get_alignment_strings) {
		if ( $get_alignment_strings / 2 % 2 ) {
			$last_attribute = "subject_alignstr";
			$xmlattr{"Hsp_hseq"} = $last_attribute;
		}
		if ( $get_alignment_strings / 4 % 2 ) {
			$last_attribute = "midline_alignstr";
			$xmlattr{"Hsp_midline"} = $last_attribute;
		}
	}

	my @hits;
	my $hit;
	my $query_id;
	my $query_definition;
	my $query_length;
	my $subject_id;
	my $subject_name;
	my $subject_definition;
	my $subject_length;

	my $linenum = 0;
	open( XML, "<$xmlfile" );

#print "$xmlfile\n";
	while ( my $line = <XML> ) {

		if ( $dbg ) { print $line }
#print $line;
		$linenum++;
		chomp $line;
		$line =~ s/[\n\r]/ /g;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		$line =~ s/&lt;/</g;
		$line =~ s/&gt;/>/g;
		$line =~ s/&amp;/&/g;
		$line =~ s/&quot;/"/g;
		$line =~ s/&apos;/'/g;

		# check XML tag
		if ( substr( $line, 0, 1 ) eq "<" ) {
			$line = substr( $line, 1 );
			my ($xml_tag) = split( ">", $line );
			my $attrname = $xmlattr{$xml_tag};
			if ( !defined $attrname ) { next }

			# get tag value
			my $attrval;
			$line = substr( $line, length($xml_tag) + 1 );
			my $close_tag = "</" . $xml_tag . ">";
			my $eod = index( $line, $close_tag );
			if ( $eod == -1 ) {
				$attrval = $line;
			}
			else {
				$attrval = substr( $line, 0, $eod );
			}

			# save tag value
			if ( $attrname eq "subject_definition" ) {
				$subject_definition = $attrval;
				$subject_definition =~ s/^[ >]*//;
				$subject_definition =~ s/ *$//;
				$$hit{subject_definition} = $subject_definition;
				$subject_id               = get_defline_id($subject_definition);
				$subject_name = get_reference_name( $subject_id );
				$$hit{subject_id}         = $subject_id;
				$$hit{subject_name}       = $subject_name;

			}
			elsif ( $attrname eq "subject_length" ) {
				$subject_length = $attrval;
				$$hit{subject_length} = $subject_length;
			}
			elsif ( $attrname eq "query_definition" ) {
				$query_definition = $attrval;
				$query_definition =~ s/^[ >]*//;
				$query_definition =~ s/ *$//;
				$$hit{query_definition} = $query_definition;
				$query_id               = get_defline_id($query_definition);
				$$hit{query_id}         = $query_id;
			}
			elsif ( $attrname eq "query_length" ) {
				$query_length = $attrval;
				$$hit{query_length} = $query_length;
			}
			else {
				$$hit{$attrname} = $attrval;
			}
			if ( !defined $attrval ) {
				my $tmpcmd = "tail +n " . ( $linenum - 5 ) . " | head -n 20 ";
#print `$tmpcmd` . "\n";
			}

			# save hit
			if ( $attrname eq $last_attribute ) {
				my $orientation = 1;
				if ( $$hit{query_frame} < 0 ) { $orientation = -1 }
				if ( $$hit{query_left} > $$hit{query_right} ) {
					$orientation = -$orientation;
					my $tmp = $$hit{query_left};
					$$hit{query_left}  = $$hit{query_right};
					$$hit{query_right} = $tmp;
				}
				if ( $$hit{subject_left} > $$hit{subject_right} ) {
					$orientation = -$orientation;
					my $tmp = $$hit{subject_left};
					$$hit{subject_left}  = $$hit{subject_right};
					$$hit{subject_right} = $tmp;
				}
				$$hit{orientation} = $orientation;
				calculate_blast_percentages( $hit );

				$$hit{has_frameshift} = 0;
#print "BLAST HSP\n";
#print_blasthits( $hit );
				# save this hit and start a new one
				push @hits, $hit;

				my %newhit;
				$newhit{query_id}           = $query_id;
				$newhit{query_definition}   = $query_definition;
				$newhit{query_length}       = $query_length;
				$newhit{subject_id}         = $subject_id;
				$newhit{subject_name}       = $subject_name;
				$newhit{subject_definition} = $subject_definition;
				$newhit{subject_length}     = $subject_length;
				$hit                        = \%newhit;
			}
		}
	}

	close XML;
	return @hits;
}

# adjust blast result coordinates back to original sequences
sub offset_adjustment {
	my ( $qstart, $qend, $qlength, $sstart, $send, $slength, @rawhsps ) = @_;

	# adjust positions for subject and query offsets
	my @hsps;
	my $qoffset = $qstart - 1;
	my $soffset = $sstart - 1;
	for my $hit (@rawhsps) {
		$$hit{query_left}    += $qoffset;
		$$hit{query_right}   += $qoffset;
		$$hit{subject_left}  += $soffset;
		$$hit{subject_right} += $soffset;

		# adjust query frame
		if ( $$hit{query_frame} > 0 ) {
			$$hit{query_frame} = ( $$hit{query_left} - 1 ) % 3 + 1;
		}
		else {
			my $reverse_start = $qlength - $$hit{query_right} + 1;
			$$hit{query_frame} = -( ( $reverse_start - 1 ) % 3 + 1 );
		}

		# adjust subject length
		if ( defined $slength ) { $$hit{subject_length} = $slength }
		calculate_blast_percentages( $hit );

		# save adjusted hit
		push @hsps, $hit;
	}

	return @hsps;
}

# combined hsps between same subject/query and in same orientation into multri-hsp "hits"
sub combine_hits {
	my ( $genome, @rawhsps ) = @_;
	my $frameshift_detection = get_parameter("frameshift_detection");
	my $max_aa_gap     = get_parameter("max_aa_gap");
	my $genome_sequence = $$genome{sequence};
	my $dbg = 0;
#	if ( @rawhsps ) {
#		if ( $rawhsps[0]{subject_id} eq "ABN10838.1" ) { $dbg = 1 }
#	}

	my %rawhits;
	for my $hsp ( @rawhsps ) {
		my $key = "$$hsp{query_id}|$$hsp{orientation}|$$hsp{subject_id}";
		if ( ! exists $rawhits{$key} ) {
			my @tmp = ( $hsp );
			$rawhits{$key} = \@tmp;
		}
		else {
			push @{$rawhits{$key}}, $hsp;
		}
	}

	my @hits;
	for my $key ( keys %rawhits ) {
		my @hsps = sort { compare_qry_positions( $a, $b ) } @{$rawhits{$key}};
		my $positions = @hsps;
		my $permutations = 1;
		for my $i ( 1..$positions ) { $permutations = $permutations * 2 }
		$permutations--;
		if ( $dbg ) {
			print "\nAVAILABLE HSPS ($permutations permutations)\n";
			print_blasthits( @hsps );
		}

		my $allow_splicing = allow_splicing( $hsps[0] );
		my $allow_slippage = allow_ribosomal_slippage( $hsps[0] );
		for my $permutation ( 1..$permutations ) {
			my $selector = $permutation;
			my @hithsps;
			for my $position ( 1..$positions ) {
				if ( $selector % 2 == 1 ) { push @hithsps, $hsps[$position-1] }
				$selector = $selector / 2;
			}

			if ( $dbg ) {
				print "permutation #$permutation\n";
				print_blasthits( @hithsps );
			}

			my $frameshifts = 0;
			my $overlaps = 0;
			for my $i ( 0..@hithsps-2 ) {
				for my $j ( $i+1..@hithsps-1 ) {
					my $qgapleft =
					  minval( $hithsps[$i]{query_right}, $hithsps[$j]{query_right} ) + 1;
					my $qgapright =
					  maxval( $hithsps[$i]{query_left}, $hithsps[$j]{query_left} ) - 1;
					my $qgapsize = $qgapright - $qgapleft + 1;

					my $sgapsize =
						maxval( $hithsps[$i]{subject_left}, $hithsps[$j]{subject_left} )
						- minval( $hithsps[$i]{subject_right}, $hithsps[$j]{subject_right} )
						- 1;

					if ( $frameshift_detection
							&& $qgapsize < 1 && $sgapsize <= 3
							&& $sgapsize >= -6 && abs( 3*$sgapsize - $qgapsize ) <= 12
							&& $hithsps[$i]{evalue} == $hithsps[$j]{evalue}
							&& $hithsps[$i]{query_frame} != $hithsps[$j]{query_frame} ) {
						if ( $dbg ) {
							print "FRAMESHIFT S $sgapsize  Q $qgapsize  E $hithsps[$i]{evalue}  F $hithsps[$i]{query_frame}|$hithsps[$j]{query_frame}\n"
						}
						$frameshifts++;
						last;
					}
					if ( -$sgapsize >= 0.5 *
								maxval( 120, minval( $hithsps[$i]{subject_coverage}, $hithsps[$j]{subject_coverage} ) )
							|| -$qgapsize >= 0.5 *
								maxval( 40, minval( $hithsps[$i]{query_coverage}, $hithsps[$j]{query_coverage} ) ) ) {
						if ( $dbg ) {
							print "OVERLAP S $sgapsize / ( $hithsps[$i]{subject_coverage}, $hithsps[$j]{subject_coverage} )  Q $qgapsize / ( $hithsps[$i]{query_coverage}, $hithsps[$j]{query_coverage} )\n"
						}
						$overlaps++;
						last;
					}
				}
				if ( $overlaps ) { last }
			}
			if ( $overlaps ) { next }
			if ( $allow_slippage && @hithsps > 2 ) { next }

			if ( ! $allow_splicing && ! $allow_slippage && @hithsps - $frameshifts > 1 ) {
				my %fusedhsp = %{$hithsps[0]};
				if ( $dbg ) {
					print "FUSE\n";
					print_blasthits( \%fusedhsp );
				}
				for my $k ( 1..@hithsps-1 ) {
					my $tmphsp = fuse_hsps( \%fusedhsp, $hithsps[$k] );
					%fusedhsp = %$tmphsp;
					if ( $dbg ) {
						print "FUSE\n";
						print_blasthits( \%fusedhsp );
					}
					#if ( $fusedhsp{query_frame} == 0 ) { last }
				}
				#$fusedhsp{has_frameshift} = 1;
				#if ( $fusedhsp{query_frame} == 0 ) { next }
				@hithsps = ( \%fusedhsp );
			}

			my %hit = %{$hithsps[0]};
			$hit{hsps} = \@hithsps;
			if ( @hithsps > 1 ) { $hit{query_frame} = 0 }
			sum_hit_hsps( \%hit );
			if ( $allow_slippage ) { $hit{has_frameshift} = 0 }
			calculate_blast_percentages( \%hit );
			push @hits, \%hit;
			if ( $dbg ) {
				print "HIT\n";
				print_blasthits( \%hit );
				print_blasthits( @{$hit{hsps}} );
			}
		}
	}

	return @hits;

}
sub join_hsps {
	my ( $genome, @rawhsps ) = @_;
	my $frameshift_detection = get_parameter("frameshift_detection");
	my $max_aa_gap     = get_parameter("max_aa_gap");
	my $genome_sequence = $$genome{sequence};
	my $dbg = 0;

	# merge overlapping hsps in same frame
	my @hsps;
	for my $i ( 0 .. @rawhsps-2 ){
		if ( ! defined $rawhsps[$i] ) { next }
		my $ihsp = $rawhsps[$i];

		for my $j ( $i+1 .. @rawhsps-1 ) {
			if ( ! defined $rawhsps[$j] ) { next }
			my $jhsp = $rawhsps[$j];

			# same pair of sequences and frame?
			if ( $$ihsp{query_id} ne $$jhsp{query_id}
					|| $$ihsp{subject_id} ne $$jhsp{subject_id}
					|| $$ihsp{query_frame} != $$jhsp{query_frame} ) {
				next;
			}

			# completely contained by better hsp?
			elsif ( $$ihsp{subject_left} >= $$jhsp{subject_left}
					&& $$ihsp{subject_right} <= $$jhsp{subject_right}
					&& $$ihsp{query_left} >= $$jhsp{query_left}
					&& $$ihsp{query_right} <= $$jhsp{query_right}
					&& $$ihsp{bit_score} > $$jhsp{bit_score} ) {
				$rawhsps[$j] = undef;
				next;
			}
			elsif ( $$jhsp{subject_left} >= $$ihsp{subject_left}
					&& $$jhsp{subject_right} <= $$ihsp{subject_right}
					&& $$jhsp{query_left} >= $$ihsp{query_left}
					&& $$jhsp{query_right} <= $$ihsp{query_right}
					&& $$jhsp{bit_score} > $$ihsp{bit_score} ) {
				$rawhsps[$i] = undef;
				last;
			}

			#overlapping?
			my $qgapleft =
			  minval( $$ihsp{query_right}, $$jhsp{query_right} ) + 1;
			my $qgapright =
			  maxval( $$ihsp{query_left}, $$jhsp{query_left} ) - 1;
			my $qgapsize = $qgapright - $qgapleft + 1;

			# merge the two hsps into one
			my $sgapsize =
				maxval( $$ihsp{subject_left}, $$jhsp{subject_left} )
				- minval( $$ihsp{subject_right}, $$jhsp{subject_right} )
				- 1;
			if ( $qgapsize <= 0 && $sgapsize <= 0 ) {
				$ihsp = merge_hsp( $genome_sequence, $ihsp, $jhsp );
				$rawhsps[$j] = undef;
				$rawhsps[$i] = $ihsp;
				next;
			}

			# no overlap, are hsps close enough to fuse?
			elsif ( $sgapsize <= $max_aa_gap && $qgapsize <= 3 * $max_aa_gap ) {

				# are there intervening stop codons?
				my $pep;
				if ( $qgapsize < 3 ) {
					$pep = "";
				}
				elsif ( $$ihsp{orientation} == -1 ) {
					$pep = DNA2AA( subsequence( $genome_sequence, $qgapright, $qgapleft ) );
				}
				else {
					$pep = DNA2AA( subsequence( $genome_sequence, $qgapleft, $qgapright ) );
				}
				if ( index( $pep, "*" ) < 0 ) {
					my $ihsp = fuse_hsps( $ihsp, $jhsp );
					$rawhsps[$j] = undef;
					$rawhsps[$i] = $ihsp;
					next;
				}
			}
		}
	}
	push @hsps,  remove_undefs( @rawhsps);
	my @tmphsps = @hsps;

	# check for cases where overlap may indicate frameshift
	for my $i ( 0 .. @tmphsps-2 ){
		my $ihsp = $tmphsps[$i];
		my $allow_slippage = allow_ribosomal_slippage( $ihsp );

		for my $j ( $i+1 .. @tmphsps-1 ) {
			my $jhsp = $tmphsps[$j];

			# same pair of sequences?
			if ( $$ihsp{query_id} ne $$jhsp{query_id}
				|| $$ihsp{subject_id} ne $$jhsp{subject_id}
				|| $$ihsp{orientation} <=> $$jhsp{orientation} ) { next }

			# calculate overlap
			my $qgapleft =
			  minval( $$ihsp{query_right}, $$jhsp{query_right} ) + 1;
			my $qgapright =
			  maxval( $$ihsp{query_left}, $$jhsp{query_left} ) - 1;
			my $qgapsize = $qgapright - $qgapleft + 1;

			my $sgapsize =
				maxval( $$ihsp{subject_left}, $$jhsp{subject_left} )
				- minval( $$ihsp{subject_right}, $$jhsp{subject_right} )
				- 1;

			# possible frameshift?
			# (blast produces two hsps in different frames with
			#  small overlap between 5' and 3' ends and same evalue)
			if ( $frameshift_detection && ! $allow_slippage
					&& $$ihsp{query_frame} != $$jhsp{query_frame}
#					&& $$ihsp{evalue} == $$jhsp{evalue}
					&& $qgapsize < 1 && $sgapsize <= 3 ) {

#				my $overlap = blasthit_span_overlap( $ihsp, $jhsp, 0 );
#				if ( $overlap > minval( $$ihsp{query_coverage}, $$jhsp{query_coverage} ) ) { next }

				if ( $$ihsp{subject_left} < $$jhsp{subject_left} ) {
					if ( $$ihsp{orientation} == 1 ) {
						if ( $$jhsp{query_left} < $$ihsp{query_left } ) { next }
					}
					else {
						if ( $$jhsp{query_right} > $$ihsp{query_right } ) { next }
					}
				}
				else {
					if ( $$ihsp{orientation} == 1 ) {
						if ( $$ihsp{query_left} < $$jhsp{query_left } ) { next }
					}
					else {
						if ( $$ihsp{query_right} > $$jhsp{query_right } ) { next }
					}
				}

				my $fusedhsp = fuse_hsps( $ihsp, $jhsp );
				$$fusedhsp{has_frameshift} = 1;
				$$fusedhsp{query_frame} = 0;
				push @hsps, $fusedhsp;
			}
		}
	}
	return @hsps;
}

# add an hsp to a hit
# (a collection of hsps sharing the same sequences in the same orientation)
sub add_hsp_to_hit {
	my ( $hit, $hsp ) = @_;

	my %orig = %$hit;

	# append to hsp lists
	push @{ $$hit{hsps} }, $hsp;

	# adjust hit statistics
	sum_hit_hsps( $hit );
	calculate_blast_percentages( $hit );

#if ( $$hit{subject_id} eq "ACN88086.1" ) {
#	print "**************************************\nADD HSP!\n";
#	print_blasthits( \%orig, $hsp, $hit, @{$$hit{hsps}} );
#	print "**************************************\n";
#}
	return $hit;
}

# create a new hsp
sub add_new_hsp {
	my ( $genome, $hit, $sbjleft, $sbjright, $qryleft, $qryright, $numid, $numsim ) = @_;

#print "add_new_hsp(" . join( ", ", ( $genome, $hit, $sbjleft, $sbjright, $qryleft, $qryright, $numid, $numsim ) ) . "\n";
#	my %orig = %$hit;
	my %new = %$hit;
	delete $new{hsps};
	delete $new{query_alignstr};
	delete $new{subject_alignstr};
	delete $new{midline_alignstr};

	$new{subject_left} = $sbjleft;
	$new{subject_right} = $sbjright;
	$new{subject_coverage} = $sbjright - $sbjleft + 1;
	$new{longest_orf} = $new{subject_coverage};
	my $ref = get_reference_seq( $$hit{subject_id} );
	$new{subject_alignstr} = substr( $$ref{sequence}, $sbjright-1, $new{subject_coverage} );

	$new{query_left} = $qryleft;
	$new{query_right} = $qryright;
	$new{query_coverage} = $qryright - $qryleft + 1;
	$new{query_alignstr} = substr( $$genome{sequence}, $qryleft-1, $new{query_coverage} );

	$new{alignment_length} = maxval( $new{subject_coverage}, $new{query_coverage} / 3 );
	$new{num_identical} = $numid;
	$new{num_similar} = $numsim;

	$new{bit_score} =
		minval(  $new{subject_coverage}, $new{query_coverage} / 3.0 ) / 4.0 + $numsim / 2.0 + $numid;
	if ( $new{bit_score} < 2 ) { $new{bit_score} = 2 }

	$new{evalue} = 5.0 - 5 * log( $new{bit_score} );

	$new{has_frameshift} = 0;

	calculate_blast_percentages( \%new );

	my @hsps = @{$$hit{hsps}};
	push @hsps, \%new;
	@hsps = sort { compare_qry_positions( $a, $b )} @hsps;
	$$hit{hsps} = \@hsps;

	sum_hit_hsps( $hit );
	calculate_blast_percentages( $hit );

#print "**************************************\nNEW HSP!\n";
#print_blasthits( \%orig, \%new, $hit );
#print "**************************************\n";
	return;
}

# calculate %identity, %coverage, etc.
sub calculate_blast_percentages {
	my ( $hsp ) = @_;

	my $frameshift = $$hsp{has_frameshift};
	if ( $frameshift ) {
		for my $attr ( "subject_coverage", "query_coverage", "num_identical", "num_similar", "longest_orf", "bit_score" ) {
			$$hsp{$attr} = int( 0.90 * $$hsp{$attr} );
		}
		$$hsp{evalue} = ( $$hsp{evalue} + 1e-70 ) / 0.90;
	}
	if ( ! defined $$hsp{subject_coverage} ) {
		$$hsp{subject_coverage}  = $$hsp{subject_right} - $$hsp{subject_left}  + 1;
	}
	if ( ! defined $$hsp{query_coverage} ) {
		$$hsp{query_coverage}  = $$hsp{query_right} - $$hsp{query_left}  + 1;
	}

	$$hsp{pct_identity} = int(
		1000. * $$hsp{num_identical} / ( $$hsp{alignment_length} ) ) / 10.;
	if ( $$hsp{pct_identity} > 100 ) { $$hsp{pct_identity} = 100 };

	$$hsp{pct_similarity} = int(
		1000. * $$hsp{num_similar} / ( $$hsp{alignment_length} ) ) / 10.;
	if ( $$hsp{pct_similarity} > 100 ) { $$hsp{pct_similarity} = 100 };

	$$hsp{pct_qcoverage} =
		int( 1000. * $$hsp{query_coverage} / $$hsp{query_length} ) / 10.;
	if ( $$hsp{pct_qcoverage} > 100 ) { $$hsp{pct_qcoverage} = 100 };
	$$hsp{pct_scoverage} =
		int( 1000. * $$hsp{subject_coverage} / $$hsp{subject_length} ) / 10.;
	if ( $$hsp{pct_scoverage} > 100 ) { $$hsp{pct_scoverage} = 100 };
	$$hsp{pct_coverage} = maxval( $$hsp{pct_qcoverage}, $$hsp{pct_scoverage} );

	$$hsp{query_gaps} =
	  	int( $$hsp{alignment_length} - $$hsp{query_coverage} / 3.0 );
	$$hsp{subject_gaps} =
		$$hsp{alignment_length} - $$hsp{subject_coverage};

	if ( defined $$hsp{query_alignstr} ) {
		my @orfs = sort { length($b) <=> length($b) } split /\*+/, $$hsp{query_alignstr};
		$$hsp{query_stops} = @orfs - 1;
		$$hsp{longest_orf} = length( $orfs[0] );
	}
	else {
		$$hsp{query_stops} = 0;
		$$hsp{longest_orf} = $$hsp{query_coverage} / 3;
	}

	$$hsp{vigor_pctsimilarity} = calculate_vigor_pctsimilarity( $hsp );
	if ( $$hsp{vigor_pctsimilarity} > 100 ) { $$hsp{vigor_pctsimilarity} = 100 };
	$$hsp{vigor_matchwgt} =
	  $$hsp{vigor_pctsimilarity} / 100.0 * $$hsp{subject_length};
}

# calculate hit statistics from constituent hsps
sub sum_hit_hsps {
	my ( $hit, @arghsps ) = @_;

	my @hsps;
	if ( @arghsps ) {
		push @hsps, @arghsps;
#print "\nSUMHITHSPS ARGHSPS\n";
	}
	else {
		push @hsps, @{$$hit{hsps}};
#print "\nSUMHITHSPS HITHSPS\n";
	}
#print_blasthits( @hsps );

	for my $hsp ( @hsps ) {
		if ( ! defined $$hsp{subject_coverage} ) {
			$$hsp{subject_coverage} = $$hsp{subject_right} - $$hsp{subject_left} + 1;
		}
		if ( ! defined $$hsp{query_coverage} ) {
			$$hsp{query_coverage} = $$hsp{query_right} - $$hsp{query_left} + 1;
		}
	}

	# hit edges
	$$hit{subject_left} = $hsps[0]{subject_left};
	$$hit{subject_right} = $hsps[0]{subject_right};
	$$hit{query_left} = $hsps[0]{query_left};
	$$hit{query_right} = $hsps[0]{query_right};
	$$hit{query_frame} = $hsps[0]{query_frame};
	$$hit{evalue} = $hsps[0]{evalue};
	$$hit{bit_score} = $hsps[0]{bit_score};

	for my $i ( 1 .. @hsps-1 ) {
		if ( ${$hsps[$i]}{subject_left} < $$hit{subject_left} ) {
			$$hit{subject_left} = $hsps[$i]{subject_left};
		}
		if ( ${$hsps[$i]}{subject_right} > $$hit{subject_right} ) {
			$$hit{subject_right} = $hsps[$i]{subject_right};
		}
		if ( ${$hsps[$i]}{query_left} < $$hit{query_left} ) {
			$$hit{query_left} = $hsps[$i]{query_left};
		}
		if ( ${$hsps[$i]}{query_right} > $$hit{query_right} ) {
			$$hit{query_right} =  $hsps[$i]{query_right};
		}
		if ( ${$hsps[$i]}{query_frame} != $$hit{query_frame} ) {
			$$hit{query_frame} = 0;
		}
		if ( ${$hsps[$i]}{evalue} < $$hit{evalue} ) {
			$$hit{evalue} =  $hsps[$i]{evalue};
		}
		if ( ${$hsps[$i]}{bit_score} > $$hit{bit_score} ) {
			$$hit{bit_score} =  $hsps[$i]{bit_score};
		}
	}

	# hit alignment counts
	my @attributes = (
		"subject_coverage", "query_coverage", "alignment_length", "num_identical",
		"num_similar", "subject_gaps", "query_gaps", "query_stops", "longest_orf",
		"has_frameshift" );
	for my $attr (@attributes) {
		$$hit{$attr} = 0;
		for my $hsp ( @hsps ) {
			$$hit{$attr} += $$hsp{$attr};
		}
	}
	if ( $$hit{subject_coverage} > $$hit{subject_right} - $$hit{subject_left} + 1 ) {
		$$hit{subject_coverage} = $$hit{subject_right} - $$hit{subject_left} + 1;
	}
	if ( $$hit{query_coverage} > $$hit{query_right} - $$hit{query_left} + 1 ) {
		$$hit{query_coverage} = $$hit{query_right} - $$hit{query_left} + 1;
	}

	for my $alignment ( "query_alignstr", "midline_alignstr", "subject_alignstr" ) {
		my $alignstr = "";
		for my $hsp ( sort { compare_qry_positions( $a, $b) } @hsps ) {
			if ( ! defined $$hsp{$alignment} ) {
				$alignstr = undef;
				last;
			}
			$alignstr .= $$hsp{$alignment};
		}
		if ( defined $alignstr ) {
			$$hit{$alignment} = $alignstr;
		}
	}

	return $hit;
}

sub calculate_vigor_pctsimilarity {
	my ( $hit ) = @_;

#	my $pct = ( $$hit{pct_identity} + $$hit{pct_similarity} ) / 2.0;

# match weight
	#my $nomatch = $$hit{subject_length} - $$hit{subject_coverage};
	#if ( $nomatch <  0 ) { $nomatch = 0 }
	my $mismatch = $$hit{subject_coverage} - $$hit{num_similar};
	my $similar = $$hit{num_similar} - $$hit{num_identical};
	my $weight = 0.10 * $mismatch + 0.60 * $similar + $$hit{num_identical}; # - $nomatch * 0.10;

# start/stop penalties
#print "$$hit{subject_id}  stops: $$hit{query_stops}  lorf: $$hit{longest_orf}  slen: $$hit{subject_length}\n";
	if ( defined $$hit{query_stops}
			&& $$hit{query_stops} > 0 ) {
		if ( ! defined $$hit{longest_orf} || $$hit{longest_orf} < 0.95 * $$hit{subject_length} ) {
#print "penalty: " . ( 0.05 * $weight ) . "\n";
			$weight = 0.90 * $weight;
		}
	}

	if ( defined $$hit{has_frameshift}
			&& $$hit{has_frameshift} > 0 ) {
		$weight = 0.90 * $weight;
	}
#print "$$hit{subject_id}  left: $$hit{subject_left}  alignstr: $$hit{query_alignstr}\n";

	if ( $$hit{subject_left} > 1
			|| ( defined $$hit{query_alignstr} && uc substr( $$hit{query_alignstr}, 0, 1 ) ne "M" ) ) {
#print "penalty: " . ( 0.05 * $weight ) . "\n";
		$weight = 0.95 * $weight;
	}
#	if ( $weight < 0.0 ) { $weight = 0.0 }

# as percentage of subject length
	my $pct = 100.0 * $weight / $$hit{subject_length};

	return $pct;
}

# fuse two very close hsps into single hsp
sub fuse_hsps {
	my ( $old, $new ) = @_;
	my %tmp    = %$old;
	my $result = \%tmp;
#	my $ori = $$old{orientation};

	# calculate size of gap separating the hsps
	my $subject_gap =
	  maxval( $$new{subject_left}, $$result{subject_left} ) -
	  minval( $$new{subject_right}, $$result{subject_right} ) - 1;
	my $query_gap =
	  ( maxval( $$new{query_left}, $$result{query_left} ) -
		  minval( $$new{query_right}, $$result{query_right} ) -
		  1 ) / 3;
	my $gap_size = maxval( $subject_gap, $query_gap );

	# sum counting metrics for the two hsps
	sum_hit_hsps( $result, $old, $new );

	# adjust for gap
	$$result{subject_gaps}     += ( $gap_size - $subject_gap );
	$$result{subject_coverage} += $subject_gap;
	$$result{query_gaps}       += ( $gap_size - $query_gap );
	$$result{query_coverage}   += $query_gap;
	$$result{alignment_length} += $gap_size;
	$$result{num_similar}      += minval( $subject_gap, $query_gap ) / 3;
	$$result{num_identical}    += minval( $subject_gap, $query_gap ) / 5;

	$$result{has_frameshift} = 0;
	if ( $$old{query_frame} ne $$new{query_frame} ) {
		$$result{has_frameshift} = 1;
		$$result{query_frame} = 0;
	}

	# recalulate percentages
	calculate_blast_percentages( $result );
#if ( $ori==-1) {
#	print "**************************************\nFUSED!\n";
#	print_blasthits( $old, $new, $result );
#	print "**************************************\n";
#}
	return $result;
}

# merge two overlapping hsps into a single hsp
sub merge_hsp {
	my ( $query_sequence, $old, $new ) = @_;
	my $vigorspace = get_parameter("vigorspace");
	my $verbose    = get_parameter("verbose");
	my $hspmerge_blastopts = get_parameter( "hspmerge_blastopts" );
#	my $ori = $$old{orientation};

	# blast subject/query regions containing the two hsps
	my $query_length = length( $query_sequence );
	my $query_left = minval( $$old{query_left}, $$new{query_left} ) - 30;
	if ( $query_left < 1 ) { $query_left = 1 }
	my $query_right = maxval( $$old{query_right}, $$new{query_right} ) + 30;
	if ( $query_right > $query_length ) { $query_right = $query_length };
	my $query_segment = substr( $query_sequence, $query_left - 1,
		$query_right - $query_left + 1 );

	open( QRY, ">$vigorspace/mrgqry" );
	print QRY ">$$old{query_definition}\n$query_segment";
	close QRY;

	my $subject_id      = $$old{subject_id};
	my $subject_length  = $$old{subject_length};
	my $subject_left    = minval( $$old{subject_left}, $$new{subject_left} ) - 10;
	if ( $subject_left < 1 ) { $subject_left = 1 }
	my $subject_right   = maxval( $$old{subject_right}, $$new{subject_right} ) + 10;
	if ( $subject_right > $subject_length ) { $subject_right = $subject_length }
	my $subject_segment = substr(
		$reference_seqs{$subject_id}{sequence},
		$subject_left - 1,
		$subject_right - $subject_left + 1
	);

	open( SUBJ, ">$vigorspace/mrgsubj" );
	print SUBJ ">$$old{subject_definition}\n$subject_segment";
	close SUBJ;

	my $cmd = "formatdb -i $vigorspace/mrgsubj -l $vigorspace/formatdb.log";
	system $cmd;

	unlink "$vigorspace/mrgxml";
	$cmd = "blastall -i $vigorspace/mrgqry -d $vigorspace/mrgsubj"
	  . " $hspmerge_blastopts -m 7 > $vigorspace/mrgxml";
	system $cmd;
#print "$cmd\n";

#if ( $ori == -1 ) {
#	print "$cmd\n";
#	print "QRY " . `cat $vigorspace/mrgqry` . "\n";
#	print "SBJ " . `cat $vigorspace/mrgsubj` . "\n";
#	print "XML " . `cat $vigorspace/mrgxml` . "\n";
#}
	system "rm $vigorspace/mrgsubj* $vigorspace/mrgqry";

	my $get_alignment_strings = 1;
	if ( defined $$old{query_alignstr} )   { $get_alignment_strings += 1 }
	if ( defined $$old{subject_alignstr} ) { $get_alignment_strings += 2 }

	# use the best of the regional blast hsps as the merged result
	my @mrghsps = sort { $$b{bit_score} <=> $$a{bit_score} }
		( parse_blastxml( "$vigorspace/mrgxml", $get_alignment_strings ) );
#if ( $ori==-1 ) {
#	print "raw merge results\n";
#	print_blasthits( @mrghsps );
#}
	@mrghsps =
		offset_adjustment( $query_left, $query_right, $$old{query_length},
			$subject_left, $subject_right, $subject_length,
			@mrghsps );
#if ( $ori==-1 ) {
#	print "adjusted merge results\n";
#	print_blasthits( @mrghsps );
#}

	my $merged = shift @mrghsps;
	$$merged{has_frameshift} = 0;

#if ( $ori==-1) {
#	print "**************************************\nMERGED!\n";
#	print_blasthits( $old, $new, $merged );
#	print "**************************************\n";
#}
	return $merged;
}

########################################
# find genomic regions likely to contain genes
sub find_candidates {
	my ( $genome ) = @_;

	my $defline  = $$genome{defline};
	my $genome_length  = length( $$genome{sequence} );

	my $tile_size = get_parameter("tile_size");
	my $tile_overlap = get_parameter("tile_overlap");
	my $candidate_blastopts = get_parameter("candidate_blastopts");
	my $candidate_db = get_parameter("candidate_refdb");
	my $reference_dbsize = get_parameter("reference_dbsize");
	my $min_candidate_pctsimilarity = get_parameter("min_candidate_pctsimilarity");
	my $min_candidate_sbjcoverage = get_parameter("min_candidate_sbjcoverage");
	my $min_candidate_uniqnucs = get_parameter("min_candidate_uniqnucs");
	my $verbose = get_parameter("verbose");

	# blast genomic sequence against reference proteins
	my @hits;
	{
		my $get_alignment_strings = 1;    # get query alignment string
		my @rawhits =
		  blast_sequence( $defline, $$genome{sequence}, $candidate_db,
			"-z $reference_dbsize $candidate_blastopts", $get_alignment_strings, $tile_size,
			$tile_overlap );

		# screen hits on minimum similarity
		for my $rawhit ( sort { $$b{pct_similarity} <=> $$a{pct_similarity} } @rawhits ) {
			if ( $$rawhit{pct_similarity} >= $min_candidate_pctsimilarity ) {
				push @hits, $rawhit;
			}
			else {
				last;
			}
		}
	}

	if ( $verbose ) {
		print "\nRAW CANDIDATE HSPS\n";
		print_blasthits( @hits );
	}

	# join multi-hsp hits into a single hit
	@hits = join_hsps( $genome, @hits );
	@hits = combine_hits( $genome, @hits );

	if ( $verbose ) {
		print "\nRAW CANDIDATE HITS\n";
		print_blasthits( @hits );
	}

	if ($verbose && @hits > 1 ) {
		print "\nCANDIDATE HITS\n";
		print_blasthits( sort { compare_qry_positions( $a, $b ) } @hits );
	}

	# resolve overlapping hits
	my $i = 0;
	while ( $i < @hits ) {
		my $ihit = $hits[$i];

		if ( ! check_coverage( $min_candidate_sbjcoverage, $ihit, $genome ) ) {
#		if ( $$ihit{pct_scoverage} < $min_candidate_sbjcoverage ) {
			splice @hits, $i, 1;
			next;
		}

		my $j = $i + 1;
		while ( $j < @hits ) {
			my $jhit = $hits[$j];

			if ( candidates_overlap( $ihit, $jhit ) ) {
				my $exception =
				  candidate_comparison_exception( $genome, $ihit, $jhit );
				if ( ! defined $exception ) {
					my $comparison =
					  default_candidate_comparison( $ihit, $jhit );
					if ( $comparison > 0 ) {
						splice @hits, $j, 1;
						next;
					}
					elsif ( $comparison < 0 ) {
						splice @hits, $i, 1;
						$i--;
						last;
					}
					else {
						splice @hits, $i, 1;
						$i--;
						last;
					}
				}
				elsif ( $exception > 0 ) {
					splice @hits, $j, 1;
					next;
				}
				elsif ( $exception < 0 ) {
					splice @hits, $i, 1;
					$i--;
					last;
				}
				else {
					$j++;
				}

			}
			else {
				$j++;
			}
		}
		$i++;
	}

	@hits = sort { order_left_to_right( $a, $b ) } @hits;

	if ( $verbose ) {
		print "\nFINAL CANDIDATES\n";
		print_blasthits( @hits );
	}
	return @hits;
}

sub default_candidate_comparison {
	my ( $ihit, $jhit ) = @_;

	if ( $$ihit{vigor_matchwgt} > $$jhit{vigor_matchwgt} ) {
		return 1;
	}
	elsif ( $$jhit{vigor_matchwgt} > $$ihit{vigor_matchwgt} ) {
		return -1;
	}
	elsif ( $$ihit{bit_score} > $$jhit{bit_score} ) {
		return 1;
	}
	elsif ( $$jhit{bit_score} > $$ihit{bit_score} ) {
		return -1;
	}
	elsif ( $$ihit{vigor_pctsimilarity} > $$jhit{vigor_pctsimilarity} ) {
		return 1;
	}
	elsif ( $$jhit{vigor_pctsimilarity} > $$ihit{vigor_pctsimilarity} ) {
		return -1;
	}
	elsif ( $$ihit{subject_id} lt $$jhit{subject_id} ) {
		return 1;
	}
	elsif ( $$jhit{subject_id} lt $$ihit{subject_id} ) {
		return -1;
	}
	else {
		return 0;
	}
}

sub default_rawgene_comparison {
	my ( $ihit, $jhit ) = @_;

	if ( 0.90 * $$ihit{vigor_pctsimilarity} > $$jhit{vigor_pctsimilarity} ) {
		return 1;
	}
	elsif ( 0.90 * $$jhit{vigor_pctsimilarity}> $$ihit{vigor_pctsimilarity}  ) {
		return -1;
	}
	elsif ( 0.90 * $$ihit{vigor_matchwgt} > $$jhit{vigor_matchwgt} ) {
		return 1;
	}
	elsif ( 0.90 * $$jhit{vigor_matchwgt} > $$ihit{vigor_matchwgt} ) {
		return -1;
	}
	if ( $$ihit{vigor_pctsimilarity} + $$ihit{vigor_matchwgt} > $$jhit{vigor_pctsimilarity} + $$jhit{vigor_matchwgt} ) {
		return 1;
	}
	elsif ( $$jhit{vigor_pctsimilarity} + $$jhit{vigor_matchwgt} > $$ihit{vigor_pctsimilarity} + $$ihit{vigor_matchwgt} ) {
		return -1;
	}
	elsif ( $$ihit{bit_score} > $$jhit{bit_score} ) {
		return 1;
	}
	elsif ( $$jhit{bit_score} > $$ihit{bit_score} ) {
		return -1;
	}
	elsif ( $$ihit{subject_id} lt $$jhit{subject_id} ) {
		return 1;
	}
	elsif ( $$jhit{subject_id} lt $$ihit{subject_id} ) {
		return -1;
	}
	else {
		return 0;
	}
}

# determine if candidates overlap sufficiently
# to be mututally exclusive
sub candidates_overlap {
	my ( $hit1, $hit2 ) = @_;

	my $ignore_strand = get_parameter("ignore_strand");
	my $overlap_threshold = get_parameter("candidate_overlap_threshold");
	my $overlap_method = get_parameter("candidate_overlap_method");
	my $overlap_ismutual = get_parameter( "candidate_overlap_ismutual" );

	# calculate span overlap
	my $dbg = 0;
	#if ( $$hit1{subject_id} eq "AAQ10563.1" || $$hit2{subject_id} eq "AAQ10563.1" ) { $dbg = 1 }

	my $overlap;
	if ( $overlap_method == 0 ) {
		$overlap = blasthit_span_overlap( $hit1, $hit2, $ignore_strand );
		if ( $dbg ) {
			print "\nOVERLAP\n";
			print_blasthits( $hit1, $hit2 );
		}
	}

	# calculate hsp overlap
	else {
		if ( $dbg ) {
			print "\nOVERLAP1\n";
			print_blasthits( @{ $$hit1{hsps} } );
			print "OVERLAP2\n";
			print_blasthits( @{ $$hit2{hsps} } );
		}
		$overlap = blasthit_hsp_overlap( $hit1, $hit2, $ignore_strand );
	}

	# compare overlap to threshold
	my $name1 = get_reference_name( $$hit1{subject_id} );
	my $name2 = get_reference_name( $$hit2{subject_id} );
	my $overlapped1 = $overlap / $$hit1{query_coverage};
	my $overlapped2 = $overlap / $$hit2{query_coverage};
	if ( $dbg ) {
		print "\nOVERLAP (method $overlap_method  mutual $overlap_ismutual  threshold $overlap_threshold)\n"
			. "  1. $name1 $overlap / $$hit1{query_coverage} = $overlapped1\n"
			. "  2. $name2 $overlap / $$hit2{query_coverage} = $overlapped2\n";
	}
	if ( $overlap <= 0 ) {
		if ( $dbg ) { print "HITS DO NOT OVERLAP\n" }
		return 0;
	}

	if ( $overlap_ismutual==2 || ( $overlap_ismutual && $name1 ne $name2 )  ) {
		if ( $overlap / $$hit1{query_coverage} < $overlap_threshold ) {
			if ( $dbg ) { print "HITS DO NOT OVERLAP\n" }
			return 0;
		}
		if ( $overlap / $$hit2{query_coverage} < $overlap_threshold ) {
			if ( $dbg ) { print "HITS DO NOT OVERLAP\n" }
			return 0;
		}
		if ( $dbg ) { print "HITS OVERLAP\n" }
		return 1;
	}
	else {
		if ( $overlap / $$hit1{query_coverage} >= $overlap_threshold ) {
			if ( $dbg ) { print "HITS OVERLAP\n" }
			return 1;
		}
		if ( $overlap / $$hit2{query_coverage} >= $overlap_threshold ) {
			if ( $dbg ) { print "HITS OVERLAP\n" }
			return 1;
		}
		if ( $dbg ) { print "HITS DO NOT OVERLAP\n" }
		return 0;
	}
}

sub blasthit_span_overlap {
	my ( $hit1, $hit2, $ignore_strand ) = @_;

	if ( $ignore_strand && $$hit1{orientation} != $$hit2{orientation} ) { return 0 }

	my $overlap = minval( $$hit2{query_right}, $$hit1{query_right} ) -
		maxval( $$hit2{query_left}, $$hit1{query_left} ) + 1;

	return $overlap;
}

sub blasthit_hsp_overlap {
	my ( $hit1, $hit2, $ignore_strand ) = @_;

	my $overlap = 0;
	for my $hsp1 ( sort { compare_qry_positions( $a, $b ) } @{ $$hit1{hsps} } )	{
		for my $hsp2 ( sort { compare_qry_positions( $a, $b ) } @{ $$hit2{hsps} } )	{
			if ( ! $ignore_strand
					&& ! $$hsp2{has_frameshift} && ! $$hsp1{has_frameshift}
					&& $$hsp2{query_frame} != $$hsp1{query_frame} ) { next }
			my $hsp_overlap = minval( $$hsp2{query_right}, $$hsp1{query_right} ) -
				maxval( $$hsp2{query_left}, $$hsp1{query_left} ) + 1;
			if ( $hsp_overlap > 0 ) { $overlap += $hsp_overlap }
		}
	}

	return $overlap;
}

# remove conflicting hsps:
#    hsps that overlap on at least one sequence
#    but not proportionally on both sequences
sub remove_conflicting_hsps {
	my ( @hsps ) = @_;

## estimate a position penalty based on weighted center of alignments
## this will be used as a tie-breaker (lowest value wins) when hsps have equal bit scoes
#	my %hits;
#	for my $hsp ( sort { $$a{subject_left} <=> $$b{subject_left} } @hsps ) {
#		my $hit_id = "$$hsp{subject_id}|$$hsp{query_id}|$$hsp{query_frame}";
#		if ( !exists $hits{$hit_id} ) {
#			$hits{$hit_id}{subject_weight} = $$hsp{bit_score};
#			$hits{$hit_id}{subject_center} =
#			  ( $$hsp{subject_left} + $$hsp{subject_right} ) / 2.0 *
#			  $$hsp{bit_score};
#			$hits{$hit_id}{query_center} =
#			  ( $$hsp{query_left} + $$hsp{query_right} ) / 2.0 *
#			  $$hsp{bit_score};
#		}
#		else {
#			$hits{$hit_id}{subject_weight} += $$hsp{bit_score};
#			$hits{$hit_id}{subject_center} +=
#			  ( $$hsp{subject_left} + $$hsp{subject_right} ) / 2.0 *
#			  $$hsp{bit_score};
#			$hits{$hit_id}{query_center} +=
#			  ( $$hsp{query_left} + $$hsp{query_right} ) / 2.0 *
#			  $$hsp{bit_score};
#		}
#	}
#	for my $hit_id ( keys %hits ) {
#		$hits{$hit_id}{subject_center} =
#		  $hits{$hit_id}{subject_center} / $hits{$hit_id}{subject_weight};
#		$hits{$hit_id}{query_center} =
#		  $hits{$hit_id}{query_center} / $hits{$hit_id}{subject_weight};
#	}
#
#	for my $hsp ( @hsps ) {
#		my $hit_id = "$$hsp{subject_id}|$$hsp{query_id}|$$hsp{orientation}";
#		my $subject_offset =
#		  ( $$hsp{subject_left} + $$hsp{subject_right} ) / 2.0 -
#		  $hits{$hit_id}{subject_center};
#		my $query_estimate =
#		  $hits{$hit_id}{query_center} + $$hsp{orientation} * $subject_offset;
#		my $position_penalty =
#		  abs( $query_estimate -
#			  ( $$hsp{query_left} + $$hsp{query_right} ) / 2.0 );
#		$position_penalty = int( 10 * $position_penalty + 0.5 ) / 10.0;
#		$$hsp{position_penalty} = $position_penalty;
#	}

	# sort hsps by sequences, bit_score, and position penalty
	@hsps = sort { conflict_compare( $a, $b ) } @hsps;

#print "==============\nREMOVE OVERLAP IN\n";
#print_blasthits( @hsps );

	#print "POS PENALTY\n";
	#print_blasthits(@hsps);
	my $i = 0;
	while ( $i < @hsps ) {
		if ( ! defined $hsps[$i] ) {
			$i++;
			next;
		}

		my $ihsp = $hsps[$i];
		my $j = $i + 1;
		while ( $j < @hsps ) {
			if ( ! defined $hsps[$j] ) {
				$j++;
				next;
			}
			my $jhsp = $hsps[$j];

			# same frame and sequence pairing?
			if (   $$ihsp{subject_id} eq $$jhsp{subject_id}
				&& $$ihsp{query_id} eq $$jhsp{query_id}
				&& $$ihsp{query_frame} == $$jhsp{query_frame} ) {

#				# keep all LARGE chunks
#				if ( $$jhsp{pct_scoverage} >= 66.7 && $$ihsp{pct_scoverage} >= 66.7 ) {
#					$j++;
#					next;
#				}

				# check overlap
				my $sleft = maxval( $$ihsp{subject_left}, $$jhsp{subject_left} );
				my $sright = minval( $$ihsp{subject_right}, $$jhsp{subject_right} );
				my $allowable_soverlap = minval( $$ihsp{subject_coverage}, $$jhsp{subject_coverage} ) / 3.0;
				my $soverlap = $sright - $sleft + 1;
				my $offset = 0;
				if ( $soverlap > 0 && $soverlap < $allowable_soverlap ) {
					if ( $$ihsp{subject_left} < $$jhsp{subject_left} ) {
						if ( $$ihsp{orientation} == 1 ) {
							$offset = $$jhsp{query_left} - $$ihsp{query_right};
						}
						else {
							$offset = $$ihsp{query_left} - $$jhsp{query_right};
						}
#print "offset2=$$ihsp{orientation} * ( $$jhsp{query_left} - $$ihsp{query_left} )=$offset\n";
					}
					else {
						if ( $$ihsp{orientation} == 1 ) {
							$offset = $$ihsp{query_left} - $$jhsp{query_right};
						}
						else {
							$offset = $$jhsp{query_left} - $$ihsp{query_right};
						}
#print "offset2=$$ihsp{orientation} * ( $$ihsp{query_left} - $$jhsp{query_left} )=$offset\n";
					}
					if ( $offset > 0 ) {
						$soverlap = 0;
					}
				}

				my $qleft    = maxval( $$ihsp{query_left}, $$jhsp{query_left} );
				my $qright   = minval( $$ihsp{query_right}, $$jhsp{query_right} );
				my $qoverlap = $qright - $qleft + 1;
				my $overlap  = maxval( $soverlap, $qoverlap );

#print "  soverlap=$soverlap  qoverlap=$qoverlap\n";
				# hsps overlap
				# drop lower scoring hsp
				if ( $soverlap > 0 || $qoverlap >= 0 ) {
#print "\nremove conflicting hsp\nqoverlap=$qoverlap  soverlap=$soverlap/$allowable_soverlap  offset=$offset\n";
#print_blasthits( $hsps[$i], $hsps[$j] );
					$hsps[$j] = undef;
				}
			}
			$j++;

		}
		$i++;
	}
	@hsps = remove_undefs( @hsps );
	return @hsps;

#	# remove hsps that are out of order
#	# group hsps by subject id
#	my %rawhits;
#	for my $hsp ( @hsps ) {
#		my $key = "$$hsp{query_id}|$$hsp{orientation}|$$hsp{subject_id}";
#		if ( exists $rawhits{$key} ) {
#			my $hithsps = $rawhits{$key};
#			push @$hithsps, $hsp;
#		}
#		else {
#			my @hithsps = ( $hsp );
#			$rawhits{$key} = \@hithsps;
#		}
#	}
#
#	# check ordering on query versus ordering against subject
#	for my $key ( sort keys %rawhits ) {
#		my $retry = 1;
##my $count =  @{$rawhits{$key}};
##print "key $key  in: $count\n";
##print_blasthits( sort { compare_qry_positions( $a, $b ) } @{$rawhits{$key}} );
#
#		while ( $retry ) {
#			my $subject_left = 0;
#			$retry = 0;
#			my $ord;
#			for my $hsp ( sort { compare_qry_positions( $a, $b ) } @{$rawhits{$key}} ) {
#
#				# if order is inconsistent, remove lowest scoring hsp and try again
#				if ( $$hsp{subject_left} <= $subject_left ) {
#					$retry = 1;
#					$ord = $hsp;
#					last;
#				}
#				$subject_left = $$hsp{subject_left};
#			}
#			if ( $retry ) {
#				my @tmp = sort { $$a{bit_score} <=> $$b{bit_score} } @{$rawhits{$key}};
#				my $del = shift @tmp;
#				$rawhits{$key} = \@tmp;
##print "\nhsp out of order, removing lowest scoring hsp\n";
##print_blasthits( $ord, $del );
#			}
##$count =  @{$rawhits{$key}};
##print " => $count";
#		}
##print "\n";
#	}
#
#	# return all remaining hsps as the final result set
#	my @final;
#	for my $key ( keys %rawhits ) {
#		if ( @{$rawhits{$key}} ) { push @final, @{$rawhits{$key}} }
#	}
#
##print "REMOVE OVERLAP OUT\n";
##print_blasthits( @hsps );
##print "==============\n";
#	return @final;
}

sub remove_undefs {
	my ( @in ) = @_;
	if ( ! @in ) { return () }

	my @out;
	for my $i ( 0 .. @in-1 ) {
		if ( defined $in[$i] ) { push @out, $in[$i] }
	}
#print "remove undefs in " . scalar @in . " out " . scalar @out . "\n";
	return @out;
}
########################################`
# examine a candidate region for potential genes
sub find_rawgenes {
	my ( $genome, $candidate ) = @_;

	my $seqlen = length( $$genome{sequence} );

	my $verbose = get_parameter( "verbose" );
#if ( $$candidate{subject_name} eq "E" ) { $verbose = 1 }
	my $vigorspace = get_parameter( "vigorspace" );
	my $small_exon_autofind = get_parameter( "small_exon_autofind" );
	my $min_pctsimilarity = get_parameter( "min_rawgene_pctsimilarity" );
	my $min_coverage = get_parameter( "min_rawgene_sbjcoverage" );
	my $rawgene_reblastopts = get_parameter( "rawgene_reblastopts" );
	my $reference_db = get_parameter( "reference_db" );
	my $overlap_threshold = get_parameter( "rawgene_overlap_threshold" );
	my $ignore_strand = get_parameter( "ignore_strand" );
	if ( ! $ignore_strand ) {
		if ( $$candidate{orientation} == -1 ) {
			$rawgene_reblastopts .= " -S 2";
		}
		else {
			$rawgene_reblastopts .= " -S 1";
		}
	}
	$rawgene_reblastopts .= " -U T";

	my $dbg = 0;
	#if ( $$candidate{subject_id} eq "AAQ10563.1" ) { $dbg = 1 }

	if ( $dbg || $verbose ) {
		print "CANDIDATES $$candidate{subject_id}\n";
		print_blasthits( $candidate, @{$$candidate{hsps}} );
	}

	# get an extended cdna sequence containing the candidate's hsps
	my $masked = mask_for_candidate( $$genome{sequence}, $candidate );
	if ( $dbg ) {
		print "masked sequence=\n$masked\n";
	}

	# re-blast extended cdna
	my @hsps = blast_sequence( $$candidate{query_id}, $masked,
				$reference_db, "$rawgene_reblastopts", 1, undef, undef, $dbg );
	if ( $dbg || $verbose ) {
		print "RAWGENE RAW HSPS\n";
		print_blasthits( @hsps );
	}

	# join the resulting hsps into hits
	@hsps = join_hsps( $genome, @hsps );
	if ( $dbg || $verbose ) {
		print "RAWGENE JOINED HSPS\n";
		print_blasthits( @hsps );
	}

	my @hits = combine_hits( $genome, @hsps );
	if ( $dbg ) {
		print "RAWGENE COMBINED HITS\n";
		print_blasthits( @hits );
	}

	if ( $verbose || $dbg ) {
		print "\nrawgene hits\n";
		print_blasthits( @hits );
	}

	# screen on coverage
	{
		my @tmp = @hits;
		@hits = ();
		for my $hit ( @tmp ) {
			if ( check_coverage( $min_coverage, $hit, $genome ) ) {
				push @hits, $hit;
			}
		}
	}

	# select hits to be used as raw genes
	my @selection;
	if ( ! @hits ) { return @selection }

	my $exception = select_rawgenes_exception( $genome, $candidate, \@hits );
	if ( defined $exception ) {
		@selection = @$exception;
	}
	else {
		my $screen = "pct_scoverage";
		my $rank = "vigor_matchwgt";
		@hits = sort { $$b{$screen} <=> $$a{$screen} } @hits;
		my $bestmatch = $min_coverage;
		if ( $hits[0]{$screen} > $bestmatch ) { $bestmatch = $hits[0]{$screen} };
		$bestmatch = 0.90 * $bestmatch;

		my @raw;
		for my $hit ( @hits ) {
			if ( check_coverage( $bestmatch, $hit, $genome ) ) {
				if ( $$hit{subject_name} eq $$candidate{subject_name} ) {
					$$hit{$rank} = 1.05 * $$hit{$rank};
				}
				push @raw, $hit;
			}
		}
		@raw = sort { $$b{$rank} <=> $$a{$rank} } @raw;
		$bestmatch = $raw[0]{$rank};
		if ( $verbose ) {
			print "screened by $screen >= $bestmatch ranked by $rank\n";
			print_blasthits( @raw );
		}

		for my $i ( 0 .. @raw-1 ) {
			if ( ! defined $raw[$i] ) { next }
			if ( $raw[$i]{$rank} < 0.90 * $bestmatch ) { last }

			for my $j ( $i+1 .. @raw-1 ) {
				if ( ! defined $raw[$j] ) { next }
				if ( $raw[$j]{$rank} < 0.90 * $bestmatch ) { last }

				my $jwidth = $raw[$j]{query_coverage};
				my $iwidth = $raw[$i]{query_coverage};
				my $overlap = blasthit_hsp_overlap( $raw[$i], $raw[$j], 0 );

				if ( $verbose ) {
					print "overlap $overlap  i=$iwidth  j=$jwidth\n";
					print_blasthits( $raw[$i], $raw[$j] );
				}
				if ( $overlap >= $overlap_threshold * $jwidth ) {
					$raw[$j] = undef;
					if ( $verbose ) { print "drop j\n" }
					next;
				}
				elsif ( $overlap >= $overlap_threshold * $iwidth
							&& 0.90 * $raw[$j]{$rank} >= $raw[$i]{$rank} ) {
					$raw[$i] = undef;
					if ( $verbose ) { print "drop i\n" }
					last;
				}
			}
			if ( defined $raw[$i] ) { push @selection, $raw[$i] }
		}
	}

	if ( $verbose || $dbg) {
		print "\n$$candidate{subject_id} SELECTED RAW GENES\n";
		print_blasthits(@selection);
	}

	# generate all possible small exon and frameshift alternatives for this rawgene
	my %finalraw;
	for my $selected ( @selection ) {
		my @tmp = ( $selected );
		{
			my $splice_permutations = find_small_exons( $genome, $selected );
			if ( defined $splice_permutations ) {
				push @tmp, @$splice_permutations;
			}
		}
		if ( $$selected{has_frameshift} && ! allow_ribosomal_slippage( $selected ) ) {
			my $frame_permutations = find_frameshift_permutations( $genome, $selected );
			if ( defined $frame_permutations ) {
				push @tmp, @$frame_permutations;
			}
		}

		# associate rawgenes sharing reference sequence to
		# prevent redundant gene calls based on lesser alternatives
		for my $tmpgene ( @tmp ) {
			my $refid = $$tmpgene{subject_id};
			if ( exists $finalraw{$refid} ) {
				if ( exists $finalraw{$refid}{permutations} ) {
					push @{$finalraw{$refid}{permutations}}, $tmpgene;
				}
				else {
					my @permutations = ( $tmpgene );
					$finalraw{$refid}{permutations} = \@permutations;
				}
			}
			else {
				$finalraw{$refid} = $tmpgene;
			}
		}
	}

	# return selected rawgenes
	my @rawgenes = sort { $$b{vigor_pctsimilarity} <=> $$a{vigor_pctsimilarity} } values %finalraw;
	return @rawgenes;
}

sub find_frameshift_permutations {
	my ( $genome, $rawgene ) = @_;

#print "MAKE FS permutationS\n";
#print_blasthits( $rawgene );

	if ( ! $$rawgene{has_frameshift} ) { return undef }
	if ( @{$$rawgene{hsps}} > 1 ) { return undef }

	my $seqlen = length( $$genome{sequence} );

	my @permutations;
	for my $shift ( -1..1 ) {
		if ( $shift != 0 ) {
			my %newgene = %$rawgene;
			$newgene{query_left} += $shift;
			if ( $newgene{query_left} < 1) { $newgene{query_left} += 3 }
			elsif ( $newgene{query_left} > $seqlen ) { $newgene{query_left} -= 3 }
			if ( exists $newgene{hsps} ) {
				my %hsp = %{${$$rawgene{hsps}}[0]};
				$hsp{query_left} = $newgene{query_left};
				my @hsps = ( \%hsp );
				$newgene{hsps} = \@hsps;
			}
			push @permutations, \%newgene;
		}
	}

#print "FS permutationS\n";
#print_blasthits( @permutations );

	return \@permutations;

}

# attempt to find small exons missed by blast
sub find_small_exons {
	my ( $genome, $hit ) = @_;

	if ( ! allow_splicing($hit) ) { return undef }
	if ( ! get_parameter("small_exon_autofind") ) { return undef }
	my $min_intron_size = get_parameter("min_intron_size");
	my $max_intron_size = get_parameter("max_intron_size");
	my $maxexpr = 5;

	my $dbg = 0;
	if ( get_parameter( "verbose" ) ) { $dbg = $$hit{subject_left} > 1 || $$hit{subject_right} < $$hit{subject_length} }
	if ( $dbg ) {
		print "find amll exons\n";
		print_blasthits( $hit );
	}

	my $donor = "(";
	for my $splice ( keys %{ get_splice_donors() } ) {
		$donor .= $splice . "|";
	}
	$donor =~ s/.$/)/;

	my $acceptor = "(";
	for my $splice ( keys %{ get_splice_acceptors() } ) {
		$acceptor .= $splice . "|";
	}
	$acceptor =~ s/.$/)/;

	my @newraws;
	my %codons = %{ get_codons() };
	#push @newraws, $hit;

	# find potential missing starts
	my @starts;
	if ( $$hit{subject_left} > 1 ) {
		my %begins;
		my $maxsz = $maxexpr;
		while ( $maxsz > 2 && ! @starts ) {
			my $subjend = $$hit{subject_left} - 2;
			if ( $subjend > $maxsz ) { $subjend = $maxsz }
			elsif ( $subjend < 2 ) { $subjend = 2 }
			if ( $subjend < $maxsz ) { $maxsz = $subjend }
			my $aaseq = uc substr( ${get_reference_seq($$hit{subject_id})}{sequence}, 0, $subjend );
			my $expression = "";
			for my $char ( split //, $aaseq ) {
				$expression .= "(";
				for my $codon ( keys %codons ) {
					if ( $codons{$codon} eq $char ) {
						$expression .= $codon . "|";
					}
				}
				$expression =~ s/.$/)/;
			}
			$expression = uc $expression;
			my $dna_start = $$hit{query_left} - $max_intron_size - 100;
			my $dna_end = $$hit{query_left} - $min_intron_size + 10;
			if ( $$hit{orientation} == -1 ) {
				$dna_start = $$hit{query_right} + $min_intron_size - 10;
				$dna_end = $$hit{query_right} + $max_intron_size + 100;
				$expression = reverse_complement( $expression );
				$expression =~ tr/[]()/][)(/;
				my $donor_string = reverse_complement( $donor );
				$donor_string =~ tr/[]()/][)(/;
				$expression = "$donor_string.{0,90}?$expression";
			}
			else {
				$expression = "$expression.{0,90}?$donor";
			}
			for my $optbase ( 1..$maxsz ) {
				my $regexp = $expression;
				if ( $optbase > 1 && $maxsz > 2 ) {
					$regexp = make_base_optional( $regexp, $optbase, $$hit{orientation} );
				}
				if ( $dbg ) { print "find start $regexp\n" }
				for my $start ( find_regexp( $regexp, $$genome{sequence}, $dna_start, $dna_end, $dbg ) ) {
					if ( $$hit{orientation} == -1 ) {
						if ( exists $begins{$$start{end}} ) {
							if ( $begins{$$start{end}}{begin} <= $$start{begin} ) { next }
						}
						$begins{$$start{end}} = $start;
					}
					else {
						if ( exists $begins{$$start{begin}} ) {
							if ( $begins{$$start{begin}}{end} >= $$start{end} ) { next }
						}
						$begins{$$start{begin}} = $start;
					}
					push @starts, $start;
				}
#				if ( $optbase==1 && @starts ) { last }
			}
			$maxsz--;
		}
		@starts = unique( @starts );
		while ( @starts > 20 ) {
			pop @starts;
		}

	}

	# find potential missing stops
	my @stops;
	if ( $$hit{subject_right} < $$hit{subject_length} ) {
		my %ends;
		my $maxsz = $maxexpr;
		while ( $maxsz > 2 ) {
			my $subjbegin = $$hit{subject_right} + 2;
			if ( $subjbegin > $$hit{subject_length}-1 ) { $subjbegin = $$hit{subject_length}-1 }
			elsif ( $subjbegin < $$hit{subject_length}-$maxsz-1 ) { $subjbegin = $$hit{subject_length}-$maxsz-1 }
			if ( $$hit{subject_length} - $subjbegin + 1 < $maxsz ) { $maxsz = $$hit{subject_length} - $subjbegin + 1 }
			my $aaseq = substr( ${get_reference_seq($$hit{subject_id})}{sequence}, $subjbegin-1 );
			$aaseq .= "*";
			my $expression = "";
			for my $char ( split //, $aaseq ) {
				$expression .= "(";
				for my $codon ( keys %codons ) {
					if ( $codons{$codon} eq $char ) {
						$expression .= $codon . "|";
					}
				}
				$expression =~ s/.$/)/;
			}
			$expression = uc $expression;
			my $dna_start = $$hit{query_right} + $min_intron_size - 10;
			my $dna_end = $$hit{query_right} + $max_intron_size + 100;
			if ( $$hit{orientation} == -1 ) {
				$dna_start = $$hit{query_left} - $max_intron_size - 100;
				$dna_end = $$hit{query_left} - $min_intron_size + 10;
				$expression = reverse_complement( $expression );
				$expression =~ tr/[]()/][)(/;
				my $acceptor_string = reverse_complement( $acceptor );
				$acceptor_string =~ tr/[]()/][)(/;
				$expression = "$expression(.){0,90}?$acceptor_string";
			}
			else {
				$expression = "$acceptor.{0,90}?$expression";
			}

			for my $optbase ( 1..$maxsz ) {
				my $regexp = $expression;
				if ( $optbase > 1 && $maxsz > 2 ) {
					$regexp = make_base_optional( $regexp, -$optbase, $$hit{orientation} );
				}
				if ( $dbg ) { print "find end for $regexp\n" }
				for my $stop ( find_regexp( $regexp, $$genome{sequence}, $dna_start, $dna_end, $dbg ) ) {
					if ( $$hit{orientation} == -1 ) {
						if ( exists $ends{$$stop{begin}} ) { next }
						$ends{$$stop{begin}} = $stop;
					}
					else {
						if ( exists $ends{$$stop{end}} ) { next }
						$ends{$$stop{end}} = $stop;
					}
					push @stops, $stop;
				}
#				if ( $optbase==1 && @stops ) { last }
			}
			$maxsz--;
		}
		@stops = unique( @stops );
		if ( @stops > 5 ) {
			splice @stops, @stops-5;
		}
	}

	for my $start ( @starts ) {
		my %tmp = %$hit;
		my $left = $$start{begin};
		my $right = $$start{begin} + 5;
		if ( $$hit{orientation} == -1 ) {
			$left = $$start{end} - 5;
			$right = $$start{end};
		}
		add_new_hsp( $genome, \%tmp, 1, 2, $left, $right, 2, 2 );
		push @newraws, \%tmp;

		if ( @stops ) {
			for my $stop ( @stops ) {
				my %tmp2 = %tmp;
				$left = $$stop{end} - 8;
				$right = $$stop{end} - 3;
				if ( $$hit{orientation} == -1 ) {
					$left = $$stop{begin} + 3;
					$right = $$stop{begin} + 8;
				}
				add_new_hsp( $genome, \%tmp2, $$hit{subject_length}-1, $$hit{subject_length}, $left, $right, 2, 2 );
#print "ADD RAWGENE\n";
#print_blasthits( @{$tmp2{hsps}} );
				push @newraws, \%tmp2;
			}
		}
	}
	for my $stop ( @stops ){
		my %tmp = %$hit;
		my $left = $$stop{end} - 8;
		my $right = $$stop{end} - 3;
		if ( $$hit{orientation} == -1 ) {
			$left = $$stop{begin} + 3;
			$right = $$stop{begin} + 8;
		}
		add_new_hsp( $genome, \%tmp, $$hit{subject_length}-1, $$hit{subject_length}, $left, $right, 2, 2 );
		push @newraws, \%tmp;
	}
	if ( $dbg ) {
		print "small exon permutations\n";
		print_blasthits( @newraws );
	}

	return \@newraws;
}

sub make_base_optional {
	my ( $regexp, $pos, $ori ) = @_;

	my $base = $ori * $pos;

	my @bases = split /\(/, $regexp;
	if ( $base > @bases ) { return $regexp }

	if ( $base > 0 ) {
		$bases[$base] = "...){0,1}";
	}
	else {
		$bases[@bases+$base+1] = "...){0,1}";
	}
	$regexp = join( "(", @bases );
	return $regexp;
}

sub find_regexp {
	my ( $regexp, $string, $start, $stop, $dbg ) = @_;

	if ( ! defined $dbg ) { $dbg = 0 }
	if ( $dbg ) { print "search for $regexp from $start-$stop\n" }

	pos($string) = maxval( $start-1, 0 );
	my @matches;
	while ( $string =~ m/$regexp/g ) {
		my %match;
		$match{string} = $&;
		$match{end} = pos($string);
		if ( $match{end} > $stop ) { last }
		$match{begin} = $match{end} - length( $match{string} ) + 1;
		if ( $dbg ) { print "$regexp matches \"$match{string}\" at $match{begin}-$match{end}\n" }
  		push @matches, \%match;
  	}

  	return @matches
}

sub mask_for_candidate {
	my ( $sequence, $candidate ) = @_;

	my $extend_span = get_parameter( "rawgene_extend_span" );
	my $extend_hsp = get_parameter( "rawgene_extend_hsp" );

	my $seqlen = length( $sequence );
	my $ori   = $$candidate{orientation};
	my $dbg = 0;
	#if ( $$candidate{subject_id} eq "AAQ10563.1" ) { $dbg = 1 }

	my $allow_splicing = allow_splicing($candidate);
	my $allow_slippage = allow_ribosomal_slippage($candidate);
	if ( $dbg ) {
		print "CANDIDATE $$candidate{subject_id}  splicing=$allow_splicing  extend=$extend_span,$extend_hsp\n";
		print_blasthits( $candidate, @{$$candidate{hsps}} );
	}

	my @hsps = sort { compare_sbj_positions( $a, $b ) } @{$$candidate{hsps}};
	if ( ! $allow_splicing && ! $allow_slippage ) {
		my %hsp = %{$hsps[0]};
		$hsp{query_left} = $$candidate{query_left};
		$hsp{query_right} = $$candidate{query_right};
		$hsp{subject_left} = $$candidate{subject_left};
		$hsp{subject_right} = $$candidate{subject_right};
		@hsps = ( \%hsp );
	}

	# initialize exons from hsps
	my @tmpexons;
	for my $i ( 0 .. @hsps-1 ) {
		my $hsp = $hsps[$i];
		my %exon;
		if ( $ori == -1 ) {
			$exon{dna_start} = $$hsp{query_right};
			$exon{dna_end} = $$hsp{query_left};
		}
		else {
			$exon{dna_start} = $$hsp{query_left};
			$exon{dna_end} = $$hsp{query_right};
		}
		if ( $dbg ) {
			print "1. ORI $ori  HSP Q: $$hsp{query_left}-$$hsp{query_right} S: $$hsp{subject_left}-$$hsp{subject_right}  EXON $exon{dna_start}-$exon{dna_end}\n";
		}

		# stretch exons to fully cover subject
		my $subject_prev;
		if ( $i > 0 ) {
			$subject_prev = ${$hsps[$i-1]}{subject_right};
		}
		else {
			$subject_prev = 0;
		}
		if ( ${$hsps[$i]}{subject_left} > $subject_prev + 1 ) {
			$exon{dna_start} -= $ori * 3 * ( ${$hsps[$i]}{subject_left} - $subject_prev - 1 );
		}

		my $subject_next;
		if ( $i < @hsps-1 ) {
			$subject_next = ${$hsps[$i+1]}{subject_left};
		}
		else {
			$subject_next = $$candidate{subject_length} + 1;
		}

		if ( ${$hsps[$i]}{subject_right} < $subject_next - 1 ) {
			$exon{dna_end} += $ori * 3 * ( $subject_next - ${$hsps[$i]}{subject_right} - 1 );
		}
		if ( $dbg ) {
			print "2. ORI $ori  HSP Q: $$hsp{query_left}-$$hsp{query_right} S: $$hsp{subject_left}-$$hsp{subject_right}  EXON $exon{dna_start}-$exon{dna_end}\n";
		}

		# extend exons
		if ( $i == 0 ) {
			$exon{dna_start} -= $ori * $extend_span;
		}
		else {
			$exon{dna_start} -= $ori * $extend_hsp;
		}
		if ( $i == @hsps-1 ) {
			$exon{dna_end} += $ori * $extend_span;
		}
		else {
			$exon{dna_end} += $ori * $extend_hsp;
		}
		if ( $dbg ) {
			print "3. ORI $ori  HSP Q: $$hsp{query_left}-$$hsp{query_right} S: $$hsp{subject_left}-$$hsp{subject_right}  EXON $exon{dna_start}-$exon{dna_end}\n";
		}

		# trim to genome
		if ( $ori == 1 ) {
			if ( $exon{dna_start} < 1 ) { $exon{dna_start} = 1 }
			if ( $exon{dna_end} > $seqlen ) { $exon{dna_end} = $seqlen }
		}
		else {
			if ( $exon{dna_end} < 1 ) { $exon{dna_end} = 1 }
			if ( $exon{dna_start} > $seqlen ) { $exon{dna_start} = $seqlen }
		}

		if ( $dbg ) {
			print "4. ORI $ori  HSP Q: $$hsp{query_left}-$$hsp{query_right} S: $$hsp{subject_left}-$$hsp{subject_right}  EXON $exon{dna_start}-$exon{dna_end}\n";
		}
		push @tmpexons, \%exon;
	}

	# finalize exons
	my @exons;
	for my $exon ( sort { $ori * $$a{dna_start} <=> $ori * $$b{dna_start} } @tmpexons ) {

		if ( ! @exons ) {
			push @exons, $exon;
		}

		# merge overlapping exons
		elsif ( $ori == 1 ) {
			if ( $$exon{dna_start} > $exons[@exons-1]{dna_end} ) {
				push @exons, $exon;
			}
			else {
				$exons[@exons-1]{dna_end} = $$exon{dna_end}
			}
		}
		else {
			if ( $$exon{dna_start} < $exons[@exons-1]{dna_end} ) {
				push @exons, $exon;
			}
			else {
				$exons[@exons-1]{dna_end} = $$exon{dna_end}
			}
		}

	}

	# use exons to construct masked genome
	my $masked = "";
	my $ptr = 0;
	for my $exon ( sort { $$a{dna_start} <=> $$b{dna_start} } @exons ) {
		if ( $dbg ) {
			print "ORI $ori  ptr $ptr  EXON $$exon{dna_start}-$$exon{dna_end}\n";
		}
		# orient to forward strand
		my $left = minval( $$exon{dna_start}, $$exon{dna_end} );
		my $right = maxval( $$exon{dna_start}, $$exon{dna_end} );

		# fill in gap between exons (masked)
		if ( $left > $ptr+1 ) {
			$masked .= lc subsequence( $sequence, $ptr+1, $left-1 );
		}
		$masked .= uc subsequence( $sequence, $left, $right );
		if ( $dbg ) { print "masked=$masked\n" }

		$ptr = $right;
	}

	# add remaining genome (masked)
	if ( $ptr < $seqlen ) {
		$masked .= lc subsequence( $sequence, $ptr+1, $seqlen );
		if ( $dbg ) { print "masked=$masked\n" }
	}
	return $masked;
}

sub gene_ref_attributes {
	return ( "ref_id", "ref_length", "ref_start", "ref_end", "ref_alignment", "pep_alignment",
				"pct_refsimilarity", "pct_refsimilarity25", "pct_refidentity",
				"pct_refcoverage", "pct_reftrunc5", "pct_reftrunc3", "match_quality",
				"gene_name", "gene_definition", "gene_quality" );
}

sub get_refstat_attributes {
	return ( "pct_refcoverage","pct_refsimilarity","pct_refsimilarity25","pct_refidentity",
				"pct_reftrunc5","pct_reftrunc3" );
}

########################################
sub find_gene {
	my ( $genome, $rawgene ) = @_;
	my $min_gene_exon_size = get_parameter("min_gene_exon_size");

	my $sequence = $$genome{sequence};
	my $seq_len = length($sequence);
	my $is_mutation = 0;
	my $warning;

	# find cds in raw gene
	my $cds = find_cds( $rawgene, $sequence );
	if ( $$cds{error} ) {
		my $gene;
		$$gene{error}   = $$cds{error};
		$$gene{warning} = $$cds{warning};
		return $gene;
	}
	$warning = $$cds{warning};

	# get protein sequence
	my $protein = $$cds{protein};

	#print "FINAL PROT =$protein\n";
	my $protein_length = length($protein);

	if ( $protein =~ /(\w\*\w)/ || !defined $protein ) {
		$is_mutation = 1;
		if ($warning) {
			$warning .= ", cdna interrupted by stop codons";
		}
		else {
			$warning = "cdna interrupted by stop codons";
		}
	}
	elsif ( $protein =~ /\*$/ ) {
		$protein_length--;
	}

	# return gene
	my %gene;
	$gene{error}       = 0;
	$gene{warning}     = $warning;
	$gene{is_mutation} = $is_mutation;

	$gene{gene_id} = "$$rawgene{query_id}.$$rawgene{gene_num}";

	$gene{dna_id}      = $$rawgene{query_id};
	$gene{orientation} = $$rawgene{orientation};

	$gene{exons}       = $$cds{exons};
	$gene{cdna}        = $$cds{cdna};
	$gene{start_codon} = $$cds{start_codon};
	$gene{start_truncated} = $$cds{start_truncated};
	$gene{stop_site}   = $$cds{stop_site};
	$gene{stop_truncated}   = $$cds{stop_truncated};
	if ( exists $$cds{translation_exception} ) {
		$gene{translation_exception} = $$cds{translation_exception};
	}
	$gene{protein}        = $protein;
	$gene{protein_length} = $protein_length;

	$gene{ref_id}          = $$rawgene{subject_id};
	$gene{ref_length}      = $$rawgene{subject_length};
	my @exons = @{$$cds{exons}};
	my $first = $exons[0];
	my $last               = $exons[@exons-1];
	$gene{ref_start}       = $$first{subject_start};
	$gene{ref_end}         = $$last{subject_end};

	$gene{gene_name}       = get_reference_name( $gene{ref_id} );
	if ( ! defined $gene{gene_name} || length( $gene{gene_name} ) > 10 ) { $gene{gene_name} = $gene{gene_id} }
	$gene{gene_definition} = get_reference_definition( $gene{ref_id} );

	return \%gene;
}

sub adjust_fragment_coordinates {
	my ( $fragment, $gene ) = @_;

	$$gene{start_codon} += $$fragment{frag_offset};
	$$gene{stop_site} += $$fragment{frag_offset};

	for my $exon ( @{$$gene{exons}} ) {
		$$exon{dna_start} += $$fragment{frag_offset};
		$$exon{dna_end} += $$fragment{frag_offset};
	}

	if ( exists $$gene{translation_exception} ) {
		$$gene{translation_exception}{begin} += $$fragment{frag_offset};
		$$gene{translation_exception}{end} += $$fragment{frag_offset};
	}
}

sub validate_gene{
	my ( $genome, $rawgene, $gene ) = @_;

	if ( ! ( $$gene{error} || $$gene{is_mutation} ) ) {
		my $ref = get_reference_seq( $$gene{ref_id} );
		my $refpep = $$ref{sequence};

		score_pep_vs_reference( $$gene{gene_id}, $$gene{protein}, $$gene{ref_id}, $refpep, $gene );
		round_refstats( $gene );

		if ( validate_reference_alignment( $genome, $gene ) ) {

			$$gene{ref_alignment} = bl2seq_alignment( $genome, $gene, get_parameter( "pep_bl2seqopts") );
			$$gene{splicing_quality} = calculate_splicing_quality( $genome, $gene );
			$$gene{gene_quality} = $$gene{match_quality} + $$gene{splicing_quality};

			# warn of frameshifts in genomic sequence
			if ( defined $rawgene && $$rawgene{has_frameshift} ) {
				my $fs_message = "possible frameshift in gene";
				if ( $$gene{warning} && index( $$gene{warning}, "%" ) < 0  ) {
					$$gene{warning} .= ", $fs_message";
				}
				else {
					$$gene{warning} = $fs_message;
				}
			}
			return;
		}
	}
	if ( $$gene{error} || $$gene{is_mutation} ) {
		my $mutation = gene2mutation( $genome, $rawgene, $$gene{warning} );
		%$gene = %$mutation;
	}
}

sub finalize_merged_partials {
	my ( $genome, $gene ) = @_;

	$$gene{ref_alignment} = bl2seq_alignment( $genome, $gene, get_parameter( "pep_bl2seqopts") );
	$$gene{splicing_quality} = calculate_splicing_quality( $genome, $gene );
	$$gene{gene_quality} = $$gene{match_quality} + $$gene{splicing_quality};
}

sub score_merged_partials {
	my ( $genome, $gene ) = @_;

	my $ref = get_reference_seq( $$gene{ref_id} );
	my $refpep = $$ref{sequence};
	my $protein = $$gene{protein};
#print "reference $$gene{ref_id}\n$refpep\nraw $$gene{gene_id} prot=$$gene{protein}\nmod prot\n$protein\n";
	score_pep_vs_reference( $$gene{gene_id}, $protein, $$gene{ref_id}, $refpep, $gene );
#print "$$gene{pep_alignment}\n";
	round_refstats( $gene );
}

# find longest cds within raw gene
sub find_cds {
	my ( $rawgene, $seq ) = @_;
	my $dbg = 0;
	#if ( $$rawgene{subject_id} eq "YP_003766.2" ) { $dbg = 1 }
	my $allow_partial_genes = get_parameter("allow_partial_genes");
	my $stop_codon_readthru = get_parameter("stop_codon_readthru");

	my $seq_len = length($seq);
	my $is_partial = 0;
	my $warning;
	my $ori = $$rawgene{orientation};
	my $ref_len = $$rawgene{subject_length};
	my $ref = get_reference_seq( $$rawgene{subject_id} );
	my %refprofile = profile_peptide( $$ref{sequence} );
	if ( $dbg ) {
		print "$$rawgene{subject_id}\n";
		print_blasthits( @{$$rawgene{hsps}} );
	}

	# extract exons from raw gene by adjusting hsps
	# to splice sites and removing internal stops
	my @exons = extract_exons_from_rawgene( $seq, $rawgene );
	if ( $dbg ) {
		my $exno = 0;
		for my $exon ( @exons ) {
			$exno++;
			if ( exists $$exon{missing_exon} ) { print "...missing...\n" }
#			else { print "$exno.  $$exon{orientation}  D: $$exon{dna_start}-$$exon{dna_end}  S: $$exon{subject_start}-$$exon{subject_end}\n" }
			else { print "$exno.  D: $$exon{dna_start}-$$exon{dna_end}  S: $$exon{subject_start}-$$exon{subject_end}\n" }
		}
	}

	# find potential start codons
	# note: start_codon specifies position of 1st base in start codon
	#
	# allow user to override start codon selection
	my $start_truncated = 0;
	my $start_position = start_codon_exception( $seq, $rawgene, \@exons );
	my @starts;
	if ( defined $start_position ) {
		my %start;
		$start{position} = $start_position;
		$start{is_truncated} = 0;
		$start{penalty} = sqrt( abs( $$rawgene{subject_left}-1 ) ) / 10.0;
		push @starts, \%start;
		if ( $dbg ) { print "0. possible start\n" }
		if ( $dbg ) { print "   position=$start{position}  truncated=$start{is_truncated}\n" }
	}
	else {
		if ( exists $exons[0]{missing_exon} ) {
			shift @exons;
			my $position = $exons[0]{dna_start};
			my %start;
			$start{position} = $position;
			$start{is_truncated} = 1;
			$start{penalty} = sqrt( abs( $$rawgene{subject_left}-1 ) ) / 100.0;
			push @starts, \%start;
			if ( $dbg ) { print "1. possible start\n" }
			if ( $dbg ) { print "   position=$start{position}  truncated=$start{is_truncated}\n" }
		}

		# look for start codon upstream
		my $startcodon_search_aarange = get_parameter( "startcodon_search_aarange" );
		my ( $uprange, $downrange ) = split /:/, $startcodon_search_aarange;
		if ( $uprange =~ /\%/ ) {
			$uprange =~ s/\%//;
			$uprange = int( $uprange / 100.0  * $ref_len + 0.5 );
		}
		if ( $uprange < 5 ) { $uprange = 5 }
		$uprange = 1 - $uprange;
		if ( $downrange =~ /\%/ ) {
			$downrange =~ s/\%//;
			$downrange = int( $downrange / 100.0  * $ref_len + 0.5 );
		}
		if ( $downrange < 5 ) { $downrange = 5 }
		$downrange = 1 + $downrange;

		my $subject_start = $exons[0]{subject_start};
		my $dna_start = $exons[0]{dna_start};
		my $aa = uc codon2aa( subsequence( $seq, $dna_start, $dna_start + 2 * $ori ) );
		while ( $subject_start >= $uprange && $aa ne "*" ) {
			if ( $aa eq "M" ) {
				if ( $dbg ) { print "upstream subj: $subject_start  min: $uprange  dna: $dna_start  aa: $aa\n" }
				my %start;
				$start{position} = $dna_start;
				$start{is_truncated} = 0;
				$start{penalty} = sqrt( abs( $subject_start-1 ) ) / 10.0;
				push @starts, \%start;
				if ( $dbg ) { print "2. possible start\n" }
				if ( $dbg ) { print "   position=$start{position}  truncated=$start{is_truncated}\n" }
			}

			my $nextstart = $dna_start - 3 * $ori;

#			if ( $dbg ) { print "next=$nextstart  ori=$ori  seqlen=$seq_len\n" }
			if ( $nextstart < 1 || $nextstart > $seq_len-2 ) {
				if ( $aa ne "M" && $allow_partial_genes ) {
					my %start;
					$start{position} = $dna_start;
					$start{is_truncated} = 1;
					$start{penalty} = sqrt( abs( $subject_start-1 ) ) / 10.0;
					push @starts, \%start;
					if ( $dbg ) { print "3. possible start\n" }
					if ( $dbg ) { print "   position=$start{position}  truncated=$start{is_truncated}\n" }
				}
				if ( $dbg ) { print "  <EOS>\n" }
				last;
			}
			my $codon = subsequence( $seq, $nextstart, $nextstart + 2 * $ori );
			$aa = uc codon2aa( $codon );
			if ( $dbg && $aa eq "*" ) { print " nextaa=$aa  codon=$codon\n" }
			$dna_start = $nextstart;
			$subject_start--;
		}

		# look for start codon downstream
		$subject_start = $exons[0]{subject_start}+1;
		$dna_start = $exons[0]{dna_start} + 3 * $ori;
		$aa = uc codon2aa( subsequence( $seq, $dna_start, $dna_start + 2 * $ori ) );
		while ( $subject_start <= $downrange && $aa ne "*" ) {
			if ( $aa eq "M" ) {
				if ( $dbg ) { print "downstream subj: $subject_start  max: $downrange  dna: $dna_start  aa: $aa\n" }
				my %start;
				$start{position} = $dna_start;
				$start{is_truncated} = 0;
				$start{penalty} = sqrt( abs( $subject_start-1 ) ) / 10.0;
				push @starts, \%start;
				if ( $dbg ) { print "4. possible start\n" }
				if ( $dbg ) { print "   position=$start{position}  truncated=$start{is_truncated}\n" }
			}

			my $nextstart = $dna_start + 3 * $ori;

#			if ( $dbg ) { print "next=$nextstart  ori=$ori  seqlen=$seq_len\n" }
			if ( $nextstart < 1 || $nextstart > $seq_len-2 ) {
				if ( $aa ne "M" && $allow_partial_genes ) {
					my %start;
					$start{position} = $dna_start;
					$start{is_truncated} = 1;
					$start{penalty} = sqrt( abs( $subject_start-1 ) ) / 10.0;
					push @starts, \%start;
					if ( $dbg ) { print "5. possible start\n" }
					if ( $dbg ) { print "   position=$start{position}  truncated=$start{is_truncated}\n" }
				}
				if ( $dbg ) { print "  <EOS>\n" }
				last;
			}
			my $codon = subsequence( $seq, $nextstart, $nextstart + 2 * $ori );
			$aa = uc codon2aa( $codon );
			if ( $dbg && $aa eq "*" ) { print " nextaa=$aa  codon=$codon\n" }
			$dna_start = $nextstart;
			$subject_start++;
		}
	}

	if ( ! @starts ) {
		if ( $dbg ) { print "NO STARTS FOUND\n" }
		my $cds;
		$$cds{error} = 1;
		if ( $$rawgene{query_stops} > 0 ) {
			$$cds{warning} = "ORF interrupted by stop codon";
		}
		elsif ( @exons > 1 || @{$$rawgene{hsps}} > 1 ) {
			$$cds{warning} = "could not find valid splicing";
		}
		else {
			$$cds{warning} = "could not find start codon";
		}
		return $cds;
	}

	# choose best start codon
	my $bestscore;
	my $beststart;
	my $bestexons;
	my $bestcdna;
	my $bestpep;
	for my $start ( @starts ) {
		my $position = $$start{position};
		my $truncated = $$start{is_truncated};

		if ( $dbg ) { print "scoring start: $position  truncated? $truncated\n" }

		my @tmpexons;
		for my $exon ( @exons ) {
			my %tmp = %$exon;
			push @tmpexons, \%tmp;
		}
		my %first = %{$tmpexons[0]};
		while ( @tmpexons > 1 ) {
			if ( ! exists $first{missing_exon} ) {
				if ( $ori == 1 && $first{dna_end} >= $position ) { last }
				elsif ( $ori == -1 && $first{dna_end} <= $position ) { last }
			}
			shift @tmpexons;
			%first = %{$tmpexons[0]};
		}
		$first{subject_start} += $ori * ( $position - $first{dna_start} ) / 3;
		$first{dna_start} = $position;
		$tmpexons[0] = \%first;

		my $stop = find_stop( $seq, $rawgene, \@tmpexons, $dbg );
		if ( $$stop{error} ) {
			if ( $dbg ) { print "   NO STOP FOUND\n"}
			next;
		}

		$$start{stop} = $stop;
		my $cdna = "";
		for my $exon ( @tmpexons ) {
			$cdna .= subsequence( $seq, $$exon{dna_start}, $$exon{dna_end} );
		}
		my $pep = DNA2AA( $cdna, $$stop{translation_exception} );
		my $score = score_profile( $pep, %refprofile );
#print "start=$position  truncated=$truncated  stop=$$stop{position}  truncated=$$stop{is_truncated}  score=$score\n";
#		my $score = quick_seqscore( $pep, $$rawgene{subject_id} );
		if ( $dbg ) {
			print "start=$position  truncated=$truncated  stop=$$stop{position}  truncated=$$stop{is_truncated}  score=$score  penalty=$$start{penalty}";
			if ( defined $$stop{translation_exception} ) { print "tr exception=$$stop{translation_exception}" }
		}
		if ( ! defined $bestscore || $score - $$start{penalty} > $bestscore ) {
			$bestscore = $score - $$start{penalty};
			$beststart = $start;
			$bestexons = \@tmpexons;
			$bestcdna = $cdna;
			$bestpep = $pep;
			if ( $dbg ) { print " ***" }
		}
		if ( $dbg ) { print "\n" }
	}

	if ( ! defined $beststart ) {
		my $cds;
		$$cds{error} = 1;
		$$cds{warning} = "could not find open reading frame";
		return $cds;
	}

	# format and return cds
	my $cds;
	$$cds{error} = 0;
	$$cds{warning} = "";
	if ($warning) { $$cds{warning} = $warning }

	$$cds{exons} = $bestexons;
	$$cds{cdna} = $bestcdna;
	$$cds{start_codon} = $$beststart{position};
	$$cds{start_truncated} = $$beststart{is_truncated};
	$$cds{stop_site}   =  $$beststart{stop}{position};
	$$cds{stop_truncated} = $$beststart{stop}{is_truncated};

	my $pep = $bestpep;
	$$cds{protein} = $pep;

	if ( defined $$beststart{stop}{translation_exception} ) {
		( $$cds{translation_exception}{begin}, $$cds{translation_exception}{end} ) =
			pep_coords_to_dna( $ori, $bestexons,
				$$beststart{stop}{translation_exception}{position},
				$$beststart{stop}{translation_exception}{position} );
		$$cds{translation_exception}{aa} =  $$beststart{stop}{translation_exception}{aa};
	}
	my $stop = index( $pep, "*" ) + 1;
	if ( $stop > 0 && $stop < length( $pep ) ) {
		$$cds{warning} = "cdna interrupted by stop codons";
		$$cds{error} = 1;
	}

	return $cds;
}

sub find_stop {
	my ( $seq, $rawgene, $exons, $dbg ) = @_;

	my $allow_partial_genes = get_parameter("allow_partial_genes");
	if ( ! defined $dbg ) { $dbg = 0 }

	# find stop codon
	# note: stop_codon specifies position of last base in last codon
	#       (including the stop codon)

	my $ori = $$rawgene{orientation};
	my $seq_len = length( $seq );
	my $ref_len = $$rawgene{subject_length};
	my $ref = get_reference_seq( $$rawgene{subject_id} );
	my %refprofile = profile_peptide( $$ref{sequence} );

	my $stop_truncated = 0;
	my $translation_exception_position;
	my $translation_exception_aa;

	# allow user to override stop codon selection
	my $stop_codon = stop_codon_exception( $seq, $rawgene, $exons );
	if ( defined $stop_codon ) {
		if ( $dbg ) { print "user defined stop codon @ $stop_codon\n" }
	}

	#  search for stop codon
	else {
		if ( exists $$exons[@$exons-1]{missing_exon} ) {
			pop @$exons;
			$stop_truncated = 1;
			if ( $dbg ) { print "last exon may be missing\n" }
		}

		# extend last exon to allow for downstream stop
		my $extension = get_parameter( "stopcodon_search_aaextension" );
		if ( $extension =~ /\%/ ) {
			$extension =~ s/\%//;
			$extension = int( $extension / 100.0  * $ref_len + 0.5 );
		}
		$extension = maxval( $extension, 25 );
		my $dna_end = $$exons[@$exons-1]{dna_end};
		my $extend_end = $dna_end + $ori * 3 * ( $ref_len - $$exons[@$exons-1]{subject_end} + $extension );
		if ( $ori == -1 && $extend_end < 1 ) {
			$stop_truncated = 1;
			while ( $extend_end < 1 ) {
				$extend_end += 3;
			}
		}
		elsif ( $ori == 1 && $extend_end > $seq_len ) {
			$stop_truncated = 1;
			while ( $extend_end > $seq_len ) {
				$extend_end -= 3;
			}
		}
		$$exons[@$exons-1]{dna_end} = $extend_end;
		if ( $dbg ) { print "stop codon searched extended to $extend_end in direction $ori\n" }

		# extract cdna from exons
		my $cdna = "";
		for my $exon ( @$exons ) {
			if ( exists $$exon{missing_exon} ) { next }

			$$exon{cdna_start} = length( $cdna ) + 1;
			my $exonseq = subsequence( $seq, $$exon{dna_start}, $$exon{dna_end} );
			$cdna .= $exonseq;
			$$exon{cdna_end} = length( $cdna );
		}

		# convert to peptide and look for stop codon
		my $pep = DNA2AA( $cdna );
		my $position = index( $pep, "*" ) + 1;
		#if ( $dbg ) { print "raw peptide=$pep\n" }
		if ( $position < 1 ) {
			if ( $stop_truncated ) {
				$stop_codon = $$exons[@$exons-1]{dna_end};
				if ( $dbg ) { print "sequence truncated before stop codon, ends at $stop_codon\n" }
			}
			else {
				if ( $dbg ) { print "no stop codon found\n" }
			}
		}

		# found stop codon
		else {
			( undef, $stop_codon ) =
				pep_coords_to_dna( $$rawgene{orientation}, $exons, $position, $position );
			$stop_truncated = 0;
			if ( $dbg ) {
				my $stop1 = $stop_codon - 2*$ori;
				my $codon = subsequence( $seq, $stop1, $stop_codon );
				print "stop @ $stop1-$stop_codon=$codon\n";
			}

			# if "read thru" enabled, look for second stop codon
			if ( allow_stopcodon_readthru( $rawgene ) ) {
				my $original_position = $position;
				my $pepscore = score_profile( substr( $pep, 0, $position-1 ), %refprofile );
#print "found stop @ $position with score $pepscore\n";
				#my $pepscore = quick_seqscore( substr( $pep, 0, $position-1 ), $$rawgene{subject_id} );
				if ( $dbg ) { print "found stop @ $position with score $pepscore\n" }
				$pep = substr( $pep, 0, $position-1 ) . "X" . substr( $pep, $position );
				$position = index( $pep, "*" ) + 1;
				if ( $position > 0 ) {

					# compare peptide with "read thru"" to original
					# keep the stop that yields the best peptide
					my $readthruscore = score_profile( substr( $pep, 0, $position ), %refprofile );
					#my $readthruscore = quick_seqscore( substr( $pep, 0, $position ), $$rawgene{subject_id} );
					if ( $dbg ) { print "readthru allows new stop @ $position with score $readthruscore" }
					if ( $readthruscore > $pepscore ) {
						$translation_exception_position = $original_position;
						$translation_exception_aa = "X";
						( undef, $stop_codon ) =
							pep_coords_to_dna( $$rawgene{orientation}, $exons, $position, $position );
						if ( $dbg ) { print " ***" }
					}
					if ( $dbg ) { print "\n" }
				}
			}
		}
	}

	# validate stop codon
	if ( ! defined $stop_codon ) {
		my $stop;
		$$stop{error} = 1;
		if ( $$rawgene{query_stops} > 0 ) {
			$$stop{warning} = "ORF interrupted by stop codon";
		}
		elsif ( @$exons > 1 || @{$$rawgene{hsps}} > 1 ) {
			$$stop{warning} = "could not find valid splicing";
		}
		else {
			$$stop{warning} = "could not find stop codon";
		}
		$$stop{warning} = "could not find stop codon";
		$$stop{is_truncated} = $stop_truncated;
		if ( $dbg ) { print "\n$$stop{warning}\n" }
		return $stop;
	}
	elsif ( $stop_truncated && ! $allow_partial_genes ) {
		my $stop;
		$$stop{error} = 1;
		$$stop{warning} = "stop codon truncated";
		$$stop{is_truncated} = $stop_truncated;
		if ( $dbg ) { print "\n$$stop{warning}\n" }
		return $stop;
	}
	else {
		my $stop;
		$$stop{error} = 0;
		$$stop{position} = $stop_codon;
		$$stop{is_truncated} = $stop_truncated;
		if ( defined $translation_exception_position ) {
			$$stop{translation_exception}{position} = $translation_exception_position;
			$$stop{translation_exception}{aa} = $translation_exception_aa;
		}
		else {
			$$stop{translation_exception} = undef;
		}
		while ( @$exons && $ori * ( $$exons[@$exons-1]{dna_end} - $stop_codon ) < 0 ) {
			pop @$exons;
		}
		if ( ! @$exons ) {
			my $stop;
			$$stop{error} = 1;
			$$stop{warning} = "could not find stop codon";
		}
		else {
			$$exons[@$exons-1]{dna_end} = $stop_codon;
		}
		if ( $dbg ) {
			print_hash( "\nFOUND STOP", $stop );
			print "\n";
		}
		return $stop;
	}
}

# derive exons from hsps
#    start with the hsp boundaries
#    removing internal stop codons
#    "polish" the boundaries to splice sites while maintaining the proper framing
sub extract_exons_from_rawgene {
	my ( $sequence, $rawgene ) = @_;

	my $splice_range = get_parameter( "splicesite_search_range" );
	my $allow_partial_genes = get_parameter("allow_partial_genes" );
	my $min_generef_pctcoverage = get_parameter("min_generef_pctcoverage");
	my $ori = $$rawgene{orientation};
	my $seqlen = length( $sequence );
	my $reflen = $$rawgene{subject_length};
	my $refseq = ${get_reference_seq($$rawgene{subject_id})}{sequence};

	my $dbg = 0;
	#if ( get_reference_name($$rawgene{subject_id}) eq "PP1ab" ) { $dbg = 1 }

	if ( $dbg ) {
		my $allow_splicing = allow_splicing($rawgene);
		print "allow_splicing=$allow_splicing\n"
			. "splice_range=$splice_range\n"
			. "allow_partial_genes=$allow_partial_genes\n"
			. "min_generef_pctcoverage=$min_generef_pctcoverage\n"
			. "ori=$ori\n"
			. "genome len=$seqlen\n"
			. "ref len=$reflen\n";
	}

	# extract "exons" resulting from ribosomal slippage
	if ( $dbg ) {
		print "has_frameshift=$$rawgene{has_frameshift}  #exons="
			. @{$$rawgene{hsps}} . "  defline=$$rawgene{subject_definition}\n";
	}
	if ( @{$$rawgene{hsps}}==2 && allow_ribosomal_slippage($rawgene) ) {
		return extract_slippage_exons( $sequence, $rawgene );
	}

	# if splicing not allowed, return a single exon
	if ( ! allow_splicing($rawgene) ) {
		my %exon;
		if ( $ori == 1 ) {
			$exon{dna_start} = $$rawgene{query_left};
			$exon{dna_end} = $$rawgene{query_right};
		}
		else {
			$exon{dna_start} = $$rawgene{query_right};
			$exon{dna_end} = $$rawgene{query_left};
		}
		$exon{cdna_start} = 1;
		$exon{cdna_end} = abs( $exon{dna_end} - $exon{dna_start} ) + 1;
		$exon{subject_start} = $$rawgene{subject_left};
		$exon{subject_end} = $$rawgene{subject_right};

		my $start = $exon{dna_start};
		while ( ( $ori == 1 && $start > 3 ) || ( $ori == -1 && $start < $seqlen - 2 ) ){
			if ( $exon{subject_start} == 1 ) { last }
			$start -= $ori * 3;
			my $aa = codon2aa( subsequence( $sequence, $start, $start + 2 * $ori ) );
			if ( $aa eq "*" ) { last }
			$exon{dna_start} = $start;
			$exon{subject_start}--;
		}

		my $end = $exon{dna_end};;
		while ( ( $ori == 1 && $end < $seqlen - 2 ) || ( $ori == -1 && $end > 3 ) ){
			if ( $exon{subject_end} == $reflen ) { last }
			$end += $ori * 3;
			my $aa = codon2aa( subsequence( $sequence, $end - 2 * $ori, $end ) );
			if ( $aa eq "*" ) { last }
			$exon{dna_end} = $end;
			$exon{subject_end}++;
		}

		my @exons = ( \%exon );
		return @exons;
	}

	if ( $dbg ) { print "EXTRACT RAWGENE HSPS\n" }
	if ( $dbg ) { print_blasthits( @{$$rawgene{hsps}} ) }

	my @hsps;
	for my $hsp ( @{$$rawgene{hsps}} ) {

		my $dna_start = $$hsp{query_left};
		my $dna_end = $$hsp{query_right};
		if ( $ori == -1 ) {
			$dna_start = $$hsp{query_right};
			$dna_end = $$hsp{query_left};
		}

		my $cdna = subsequence( $sequence, $dna_start, $dna_end );
		my $aa = DNA2AA( $cdna );
		if ( $dbg ) {
			print "$dna_start-$dna_end = $aa\n";
		}
		my @frags = split /(\*+)/, $aa;
		if ( @frags == 1 ) {
			push @hsps, $hsp;
			if ( $dbg ) { print "hsp # " . scalar @hsps . " S: $$hsp{subject_left} - $$hsp{subject_right}  Q: $$hsp{query_left} - $$hsp{query_right}\n" }
		}
		else {
			my $ratio = ( $$hsp{subject_right} - $$hsp{subject_left} + 1.0 )
				/ ( $$hsp{query_right} - $$hsp{query_left} + 1.0 );
			my $dna_base = $dna_start;
			my $subject_base = $$hsp{subject_left};
			for my $frag ( @frags ) {
				if ( $frag =~ /\*/ ) {
					$dna_start += $ori * 3 * length( $frag );
				}
				else {
					$dna_end = $dna_start + $ori * 3 * length( $frag ) - $ori;
					my $subject_left = int( $subject_base + $ratio * abs( $dna_start - $dna_base ) + 0.5 );
					my $subject_right = int( $subject_base + $ratio * abs( $dna_end - $dna_base ) + 0.5 );
					my %new = %$hsp;
					$new{subject_left} = $subject_left;
					$new{subject_right} = $subject_right;
					$new{query_left} = minval( $dna_start, $dna_end );
					$new{query_right} = maxval( $dna_start, $dna_end );
					if ( @hsps > 1 ) { pop @hsps }
					push @hsps, \%new;
					if ( $dbg ) { print "hsp # " . scalar @hsps . " S: $$hsp{subject_left} - $$hsp{subject_right}  Q: $$hsp{query_left} - $$hsp{query_right}\n" }
					$dna_start += $ori * 3 * length( $frag );
				}
			}
		}

		if ( $dbg ) { print "hsp # " . scalar @hsps . " S: $$hsp{subject_left} - $$hsp{subject_right}  Q: $$hsp{query_left} - $$hsp{query_right}\n" }
	}
	my $hsp = pop @hsps;
	if ( $ori == -1 ) {
		$$hsp{query_left} = $$rawgene{query_left};
	}
	else {
		$$hsp{query_right} = $$rawgene{query_right};
	}
	push @hsps, $hsp;

	# reformat hsps as exons
	# change left/right to start/end
	@hsps = sort { compare_qry_positions( $a, $b ) } @hsps;
	my @exons;
	if ( $dbg ) { print "\nRAWEXONS\n" }
	for my $hsp (@hsps) {

		my %exon;

		if ( $ori == -1 ) {
			$exon{dna_start} = $$hsp{query_right};
			$exon{dna_end}   = $$hsp{query_left};
		}
		else {
			$exon{dna_start} = $$hsp{query_left};
			$exon{dna_end}   = $$hsp{query_right};
		}

		$exon{subject_start} = $$hsp{subject_left};
		$exon{subject_end} = $$hsp{subject_right};

		push @exons, \%exon;
		if ( $dbg ) { print "  Q: $exon{dna_start}-$exon{dna_end}  S: $exon{subject_start}-$exon{subject_end}\n" }
	}

	# are we missing the first exon?
	if ( $allow_partial_genes ) {
		my $first = $exons[0];
		my $subject_minstart = $$first{subject_start};
		my $genome_start = 1;
		my $last = $exons[ @exons-1 ];
		my $subject_maxend = $$last{subject_end};
		my $genome_end = $seqlen;
		if ( $ori == -1 ) {
			$genome_start = $seqlen;
			$genome_end = 1;
		}

		if ( $subject_minstart > 4 && $ori * ( $$first{dna_start} - $genome_start ) > 200 ) {
			my %missing;
			$missing{missing_exon} = 1;
			$missing{dna_end} = $genome_start - $ori;
			$missing{dna_start} = $missing{dna_end} + $ori;
			$missing{subject_start} = 1;
			$missing{subject_end} = $subject_minstart -1 ;
			unshift @exons, \%missing;

			if ( $dbg ) {
				print "  Missing First Exon Genome Starts at $genome_start "
				. "  Known subject Q: $$first{dna_start}-$$last{dna_end}  S: $subject_minstart-$subject_maxend"
				. "  Filler exon   Q: $missing{dna_start}-$missing{dna_end}  XS: $missing{subject_start}-$missing{subject_end}\n";
			}
		}

		# are we missing last exon
		if ( $subject_maxend < $reflen-3 && $ori * ( $genome_end - $$last{dna_end} ) > 200 ) {
			my %missing;
			$missing{missing_exon} = 1;
			$missing{dna_start} = $genome_end + $ori;
			$missing{dna_end} = $missing{dna_start} - $ori;
			$missing{subject_start} = $subject_maxend + 1;
			$missing{subject_end} = $reflen;
			push @exons, \%missing;

			if ( $dbg ) {
				print "  Missing Last Exon Genome Ends @ $genome_end"
				. "  Known subject Q: $$first{dna_start}-$$last{dna_end}  S: $subject_minstart-$subject_maxend "
				. "  filler exon   Q: $missing{dna_start}-$missing{dna_end}  S: $missing{subject_start}-$missing{subject_end}\n";
			}
		}
	}

	# expand exons to provide some flexibility for splicing
	# do not introduce a stop codon
	# do not extend beyond (theoretical) subject boundaries
	my $expansion_overlap = 10;
	for my $i ( 0 .. @exons-1 ) {
		my $exon = $exons[$i];
		if ( exists $$exon{missing_exon} ) { next }
		if ( $dbg ) { print "\n S: $$exon{subject_start}-$$exon{subject_end}  D: $$exon{dna_start}-$$exon{dna_end}\n" }
		$$exon{expanded_dna_start} = $$exon{dna_start};
		$$exon{expanded_dna_end} = $$exon{dna_end};
		$$exon{expanded_subject_start} = $$exon{subject_start};
		$$exon{expanded_subject_end} = $$exon{subject_end};
		if ( $i > 0 ) {
			my $preexon = $exons[$i-1];
			if ( ! exists $$preexon{missing_exon}
					&& $$preexon{subject_end}-$expansion_overlap < $$exon{subject_start} ) {
				my $adjust_subject = $$preexon{subject_end} - $expansion_overlap - $$exon{subject_start};
				my $adjust_dna = $ori * 3 * $adjust_subject;
				$$exon{expanded_subject_start}+= $adjust_subject;
				$$exon{expanded_dna_start}+= $adjust_dna;
			}
		}
		if ( $i < @exons-1 ) {
			my $postexon = $exons[$i+1];
			if ( ! exists $$postexon{missing_exon}
					&& $$postexon{subject_start}+$expansion_overlap > $$exon{subject_end} ) {
				my $adjust_subject = $$postexon{subject_start}+$expansion_overlap - $$exon{subject_end};
				my $adjust_dna = $ori * 3 * $adjust_subject;
				$$exon{expanded_subject_end}+= $adjust_subject;
				$$exon{expanded_dna_end}+= $adjust_dna;
			}
		}

		my %bkup = %$exon;
		$$exon{bkup} = \%bkup;
	}

	# adjust exon boundaries to splice sites
	# for each intron
	my $ref = get_reference_seq( $$rawgene{subject_id} );
	my %refprofile = profile_peptide( $$ref{sequence} );
	my $min_exon_size   = maxval( get_parameter("min_exon_size"),   1 );
	my $min_intron_size = maxval( get_parameter("min_intron_size"), 4 );
	my $exno = 1;
	while ( $exno < @exons ) {
		my $exon5 = $exons[ $exno-1 ];
		if ( exists $$exon5{missing_exon} ) {
			$exno++;
			next;
		}
		my $exon3 = $exons[ $exno ];
		if ( exists $$exon3{missing_exon} ) {
			$exno++;
			next;
		}
		if ( $dbg ) { print "$exno. find intron between $$exon5{dna_start}-$$exon5{dna_end} and $$exon3{dna_start}-$$exon3{dna_end}\n" }
		# establish potential range for intron
		my $start = $$exon5{dna_start} + $ori * $min_exon_size;
		my $end   = $$exon3{dna_end} - $ori * $min_exon_size - $ori;
		if ( $dbg ) { print "   intron between $start-$end\n" }
		# not enough room for intron, dropped smaller exon
		if ( $ori * ( $end - $start ) + 1 < $min_intron_size ) {
			if (
				abs( $$exon5{dna_end} - $$exon5{dna_start} ) >
				abs( $$exon3{dna_end} - $$exon3{dna_start} ) )
			{
				splice @exons, $exno, 1;
				next;
			}
			else {
				splice @exons, $exno-1, 1;
				next;
			}
		}

		# find all possible donor and acceptor sites
		my %donorlist = %{get_splice_donors()};
		my %acceptorlist = %{get_splice_acceptors()};
		my @donors;
		my @acceptors;
		my $donor_frame_start  = $$exon5{dna_start};
		my $acceptor_frame_end = $$exon3{dna_end};
		my $site = $start;
		while ( $ori * ( $end - $site ) > 0  ) {
			my $seq = uc subsequence( $sequence, $site, $site + $ori );
			if ( exists $donorlist{$seq}
					&& minval( abs( $site - $$exon5{dna_end} ),
							   abs( $site - $$exon5{expanded_dna_end} ) ) <= $splice_range ) {
				my %splice_site;
				$splice_site{seq} = $seq;
				$splice_site{site} = $site;
				$splice_site{framing} = abs( $site - $donor_frame_start ) % 3;
				push @donors, \%splice_site;
				if ( $dbg ) { print "D: $seq $site-" . ( $site + $ori ) . "\n" }
			}
			if ( exists $acceptorlist{$seq}
					&& minval( abs( $$exon3{dna_start} - $site ),
					           abs( $$exon3{expanded_dna_start} - $site ) ) <= $splice_range ) {
				my %splice_site;
				$splice_site{seq} = $seq;
				$splice_site{site} = $site + $ori;
				$splice_site{framing} = abs( $site + $ori - $acceptor_frame_end ) % 3;
				push @acceptors, \%splice_site;
				if ( $dbg ) { print "A: $seq $site-" . ( $site + $ori ) . "\n" }
			}
			$site += $ori;
		}

		# if the exons are in same frame with no intervening
		# stop codon we may not need an intron
		my $bestdonor;
		my $bestacceptor;
		my $bestscore;
		if ( abs( $$exon5{dna_start} - $$exon3{dna_start} ) % 3 == 0 ) {
			my $dna = subsequence( $sequence, $$exon5{dna_start}, $$exon3{dna_end} );
			my $aa = DNA2AA( $dna );
			if ( $dbg ) { print "aa=$aa\n" }
			if ( index( $aa, "*" ) < 0 ) {

				# score sequence resulting from fusing the exons
				if ( $dbg ) { print "$$exon5{dna_start}-$$exon3{dna_end}: " }
				$bestscore = score_profile( $aa, %refprofile );
				#$bestscore = quick_seqscore( $aa, $$rawgene{subject_id} );
				if ( $dbg ) { print "score $$exon5{dna_start}-$$exon3{dna_end} = $bestscore\n" }
				if ( $dbg ) { print "   no intron score: $bestscore\n" }
			}
		}

		# find best donor/acceptor combination
		for my $donor (@donors) {
		if ( $dbg ) { print "   donor $$donor{site} framing $$donor{framing}\n" }

			for my $acceptor (@acceptors) {
				if ( $dbg ) { print "     acceptor $$acceptor{site} framing $$acceptor{framing}" }

				# are donor and acceptor properly paired?
				if ( ! exists $splice_pairs{$$donor{seq}}{$$acceptor{seq}} ) {
				if ( $dbg ) { print " cannot pair $$donor{seq}-$$acceptor{seq}\n" }
					next;
				}

				# are donor/acceptor far enough apart to be intron?
				my $intron_sz = $ori * ( $$acceptor{site} - $$donor{site} ) + 1;
				if ( $dbg ) { print "$$donor{site}-$$acceptor{site} size = $intron_sz vs min $min_intron_size" }
				if ( $intron_sz < $min_intron_size ) {
				if ( $dbg ) { print " too close\n" }
					next;
				}
				if ( $dbg ) { print "\n" }
				# do donor/acceptor preserve translation frame?
				if ( ( $$donor{framing} + $$acceptor{framing} ) % 3 != 0 ) {
				if ( $dbg ) { print " framing error\n" }
					next;
				}

				# does splicing introduce a stop codon?
				my %tmp5 = %$exon5;
				$tmp5{dna_end} = $$donor{site} - $ori;
				my %tmp3 = %$exon3;
				$tmp3{dna_start} = $$acceptor{site} + $ori;

				my $aa = DNA2AA(
					subsequence( $sequence, $tmp5{dna_start}, $tmp5{dna_end} )
					  . subsequence( $sequence, $tmp3{dna_start}, $tmp3{dna_end} )
				);

				if ( $exno == @exons-1 ) { $aa =~ s/\*$// }
				if ( index( $aa, "*" ) >= 0 ) {
				if ( $dbg ) { print " introduces stop codon\n" }
					next;
				}

				# score the resuting sequence
				my $score = score_profile( $aa, %refprofile );
				#my $score = quick_seqscore( $aa, $$rawgene{subject_id} );
				$score = $score	- 100.0 * $splice_pairs{$$donor{seq}}{$$acceptor{seq}} / $reflen;

				# keep best scoring splice
				if ( ! defined $bestscore || $score > $bestscore ) {
				if ( $dbg ) { print "*** " }
					$bestscore = $score;
					$bestdonor = $donor;
					$bestacceptor = $acceptor;
				}
				if ( $dbg ) { print "$tmp5{dna_start}-$tmp5{dna_end} + $tmp3{dna_start}-$tmp3{dna_end}: $score\n" }
				if ( $dbg ) { print "gene=$aa\n ref=$refseq\n" }
			}
		}

		# no acceptable splicing, dropped smaller exon
		if ( ! defined $bestscore || $bestscore == 0 ) {
			if (
				abs( $$exon5{dna_end} - $$exon5{dna_start} ) >
				abs( $$exon3{dna_end} - $$exon3{dna_start} ) ) {
				if ( $dbg ) { print "   drop exon3\n" }

				# delete exon3
				splice @exons, $exno, 1;
				next;
			}
			else {
				if ( $dbg ) { print "   drop exon5\n" }

				# delete exon5
				splice @exons, $exno-1, 1;

				# if exon5 was not the first exon we need to restore and re-splice
				if ( $exno > 1 ) {
					$exno--;
					$exon5 = $exons[$exno];
					my %tmp = %{$$exon5{bkup}};
					$exon5 = \%tmp;
					my %bkup = %$exon5;
					$$exon5{bkup} = \%bkup;

					$exons[$exno-1] = $exon5;
				}
			}
		}

		# best alternative is NO intron, fuse the exons
		elsif ( !defined $bestdonor ) {
			if ( $dbg ) { print "   join exons\n" }
			$$exon5{dna_end} = $$exon3{dna_end};
			$$exon5{subject_end} = $$exon3{subject_end};
			splice @exons, $exno, 1;
			next;
		}

		# modify exons using the best splicing sites
		else {
			if ( $dbg ) { print "   using intron $$bestdonor{site}-$$bestacceptor{site} to modify\n" }
			if ( $dbg ) { print "   X5 Q: $$exon5{dna_start}-$$exon5{dna_end}  S: $$exon5{subject_start}-$$exon5{subject_end}\n" }
			if ( $dbg ) { print "   X3 Q: $$exon3{dna_start}-$$exon3{dna_end}  S: $$exon3{subject_start}-$$exon3{subject_end}\n" }

			# alter exon coordinates to accomodate splicing sites
			my $new_start = $$bestacceptor{site} + $ori;
			$$exon3{subject_start} += int ( $ori * ( $new_start - $$exon3{dna_start} ) / 3.0 + 0.5 );
			$$exon3{dna_start} = $new_start;

			my $new_end = $$bestdonor{site} - $ori;
			$$exon5{subject_end} += int ( $ori * ( $new_end - $$exon5{dna_end} ) / 3.0 + 0.5 );
			$$exon5{dna_end}   = $new_end;
			$exno++;
			if ( $dbg ) { print "yielding\n   X5 Q: $$exon5{dna_start}-$$exon5{dna_end}  S: $$exon5{subject_start}-$$exon5{subject_end}\n" }
			if ( $dbg ) { print "   X3 Q: $$exon3{dna_start}-$$exon3{dna_end}  S: $$exon3{subject_start}-$$exon3{subject_end}\n" }
			if ( $dbg ) { print "---------------\n" }
		}
	}

	# calculate exon cdna positions
		if ( $dbg ) { print "\nFINAL RAW EXONS--------------\n" }
		my $cdna_length = 0;
		for my $exon ( @exons ) {
		if ( exists $$exon{missing_exon} ) {
			if ( $dbg ) {
				print "  DNA: $$exon{dna_start}-$$exon{dna_end}  CDNA: missing"
					. "  SUBJ: $$exon{subject_start}-$$exon{subject_end}\n";
			}
		}
		else {
			$$exon{cdna_start} = $cdna_length + 1;
			$cdna_length += abs( $$exon{dna_end} - $$exon{dna_start} ) + 1;
			$$exon{cdna_end} = $cdna_length;
			delete $$exon{bkup};

# remove expansion attributes
			delete $$exon{expanded_dna_start};
			delete $$exon{expanded_dna_end};
			delete $$exon{expanded_subject_start};
			delete $$exon{expanded_subject_end};

			if ( $dbg ) {
				print "  DNA: $$exon{dna_start}-$$exon{dna_end}  CDNA: $$exon{cdna_start}-$$exon{cdna_end}"
					. "  SUBJ: $$exon{subject_start}-$$exon{subject_end}\n";
			}
		}
	}
	if ( $dbg ) { print "\n---------------------------------\n" }

	# return results
	return @exons;
}

# return two "ribosomal slippage exons," one pre-frameshift, one post-frameshift
sub extract_slippage_exons {
	my ( $sequence, $rawgene ) = @_;

	my $dbg = 0;
	#if ( get_reference_name($$rawgene{subject_id}) eq "PP1ab" ) { $dbg = 1 }

	# get reference sequence for comparison
	my $ref = get_reference_seq( $$rawgene{subject_id} );
	my %refprofile = profile_peptide( $$ref{sequence} );

	# find range where slippage could have occurred
	my @hsps = sort { compare_qry_positions( $a, $b ) } @{$$rawgene{hsps}};
	my $pep;
	my $left;
	my $right;
	my $offset;
	if ( $$rawgene{orientation} == 1 ) {
		$pep = DNA2AA( subsequence( $sequence, $$rawgene{query_left}, $$rawgene{query_right} ) );
		$left = $hsps[1]{query_left}-12;
		$right = $hsps[0]{query_right}+12;
		$offset = $hsps[0]{query_left};
		if ( $hsps[1]{subject_left} > $hsps[0] {subject_right}+1 ) {
			$left -= 3 * ( $hsps[0] {subject_right} - $hsps[1]{subject_left} - 1);
			$right += 3 * ( $hsps[0] {subject_right} - $hsps[1]{subject_left} - 1);
		}
	}
	else {
		$pep = DNA2AA( subsequence( $sequence, $$rawgene{query_right}, $$rawgene{query_left} ) );
		$left = $hsps[0]{query_left}-12;
		$right = $hsps[1]{query_right}+12;
		$offset = $hsps[0]{query_right};
		if ( $hsps[1]{subject_left} > $hsps[0] {subject_right} ) {
			$left -= 3 * ( $hsps[0] {subject_right} - $hsps[1]{subject_left} - 1);
			$right += 3 * ( $hsps[0] {subject_right} - $hsps[1]{subject_left} - 1);
		}
	}
	if ( $left > $right ) {
		( $right, $left ) = ( $left-6, $right+6 );
	}

	# dont look past orf A stop codon
	my @tmp = split /\*/, $pep;
	while ( @tmp ) {
		$pep = shift @tmp;
		$offset += $$rawgene{orientation} * 3 * ( length( $pep ) + 1 );
#print "L-R $left-$right  offset=$offset\n";
		if ( $offset >= $left - 100 && $offset <= $right + 100 ) { last }
	}
	if ( $$rawgene{orientation} == 1 ) {
		if ( $right >= $offset ) { $right = $offset - 1 }
	}
	else {
		if ( $left <= $offset ) { $left = $offset + 1 }
	}
	if ( $left < 1 ) { $left = 1 }
	if ( $right > length( $sequence ) ) { $right = length( $sequence ) }

	# try the possible points of slippage,
	# keep the one that yields the closest match to the reference
	my $bestscore;
	my @exons;
	for my $slip ( $left..$right ) {
		my $subject_mid = $hsps[0]{subject_right};
		if ( $$rawgene{orientation} == 1 ) {
			$subject_mid += int( ( $slip - $hsps[0]{query_right} ) / 3.0 );
		}
		else {
			$subject_mid += int( ( $hsps[0]{query_left} - $slip ) / 3.0 );
		}
		for my $fs( 0..2 ) {
			my $exon1;
			my $exon2;
			$$exon1{subject_start} = $$rawgene{subject_left};
			$$exon1{subject_end} = $subject_mid;
			$$exon2{subject_start} = $subject_mid;
			$$exon2{subject_end} = $$rawgene{subject_right};
			if ( $$rawgene{orientation} == 1 ) {
				$$exon1{dna_start} = $hsps[0]{query_left};
				$$exon1{dna_end} = $slip + $fs;
				$$exon2{dna_start} = $slip;
				$$exon2{dna_end} = $hsps[1]{query_right};
			}
			else {
				$$exon1{dna_start} = $hsps[0]{query_right};
				$$exon1{dna_end} = $slip - $fs;
				$$exon2{dna_start} = $slip;
				$$exon2{dna_end} = $hsps[1]{query_left};
			}
			my $framing = abs( $$exon1{dna_end} - $$exon1{dna_start} ) + 1
				+ abs( $$exon2{dna_end} - $$exon2{dna_start} ) + 1;
			$framing = $framing % 3;
			if ( $dbg ) {
				print "slip=$slip  fs=$fs  framing=$framing  exons=$$exon1{dna_start}..$$exon1{dna_end},$$exon2{dna_start}..$$exon2{dna_end}";
			}
			if ( ( 2 + abs( $$exon1{dna_end} - $$exon1{dna_start} + $$exon2{dna_end} - $$exon2{dna_start} ) )
					% 3 != 0 ) {
				if ( $dbg ) { print "  out-of-frame\n" }
				next;
			}
			my $pep = DNA2AA(
						subsequence( $sequence, $$exon1{dna_start}, $$exon1{dna_end} )
						. subsequence( $sequence, $$exon2{dna_start}, $$exon2{dna_end} ) );
			my $stop = index( $pep, "*" );
			if (  $stop >= 0 ) {
				if ( $stop < 0.50 * $$rawgene{subject_length} ) {
					if ( $dbg ) { print "  interrupted by stop\n" }
					next;
				}
				else {
					$pep = substr( $pep, 0, $stop );
				}
			}

			my $pepscore = score_profile( $pep, %refprofile );
			#my $pepscore = quick_seqscore( $pep, $$rawgene{subject_id} );
			if ( $dbg ) { print "  pepscore=$pepscore\n" }
			if ( ! defined $bestscore || $pepscore > $bestscore ) {
				$bestscore = $pepscore;
				@exons = ();
				push @exons, $exon1;
				push @exons, $exon2;
			}
		}
	}

	# format and return exons
	my $cdna_length = 0;
	for my $exon ( @exons ) {
		$$exon{cdna_start} = $cdna_length + 1;
		$cdna_length += abs( $$exon{dna_end} - $$exon{dna_start} ) + 1;
		$$exon{cdna_end} = $cdna_length;
	}
	return @exons;
}

# quick and dirty rotines for comparing two proteins based on kmer counts
# generate a protein kmer profile
sub profile_peptide {
	my ( $peptide ) = @_;
	my $binsz = 3;

	my $peplen = length( $peptide );
	if ( $peptide =~ /\*$/ ) { $peplen-- }

	my $pep = uc substr( $peptide, minval( 5, $peplen ) );
#	if ( substr( $pep, 0, 1 ) eq "M" ) { $pep =~ s/^M/mM/ }
	my %profile;

	$profile{aa_len} = $peplen;

#	my $vigorspace = get_parameter( "vigorspace" );
#	$profileid++;
#	$profile{file} = "$vigorspace/$profileid.fasta";
#	unlink $profile{file};
#	open( PROFILE, ">$profile{file}" );
#	print PROFILE ">$profileid\n$ucpep\n";
#	close PROFILE;

	for my $sz ( 1..$binsz ) {
		for my $i ( 0..length( $pep ) - $sz ) {
			my $bin = substr( $pep, $i, $sz );
			if ( exists $profile{bins}{$bin} ) {
				$profile{bins}{$bin}++;
			}
			else {
				$profile{bins}{$bin} = 1;
			}
		}
	}

	return %profile;
}

# compare a protein to a previous computed kmer profile
sub score_profile {
	my ( $peptide, %refprofile ) = @_;
	my $dbg = $refprofile{dbg};
	if ( ! defined $dbg ) { $dbg = 0 }
	my $binsz = 3;

	if ( $dbg ) { print "peptide=$peptide\n" }
	my %profile = profile_peptide( $peptide );

#	my $profdir = "/home/gsims/usr/bin";
#	my $cmd = "$profdir/ffpaa -l 3 -d $refprofile{file} $profile{file} | $profdir/ffpcol -a |  $profdir/ffpjsd --dice -s -r 1";
#	my $rawscore = `$cmd`;
#	my @tmp = split /  */, $rawscore;
#	shift @tmp;
#
#	my $score = $tmp[0] + ( 1.0 - abs( 1.0 - $profile{protlen}/$refprofile{protlen} ) );
#	$score = 100.0 * more FluC.cd$score / 2.0;

#print "profcmd=$cmd\nprofraw=$rawscore\nprofscore=$score\n";

	my @bins = unique( keys %{$profile{bins}}, keys %{$refprofile{bins}} );

	my $match = 0;
	my $mismatch = 0;
	for my $bin ( @bins ) {
		if ( $bin eq "aa_len") { next }
		my $pep = $profile{bins}{$bin};
		if ( ! defined $pep ) { $pep = 0 }
		my $ref = $refprofile{bins}{$bin};
		if ( ! defined $ref ) { $ref = 0 }
		$match += sqrt( length( $bin ) ) * minval( $pep, $ref );
		$mismatch += abs( $pep - $ref ) / sqrt( length( $bin ) );
	}
#	my $simscore = $match / ( $match + $mismatch );
#	my $lenscore = 1.0 - abs( 1.0 - $profile{aa_len} / $refprofile{aa_len} );
#	return 100.0 * $lenscore * $simscore;
	my $len_mismatch = abs( $profile{aa_len} / $refprofile{aa_len} );
	my $score = 100.0 * $match / ( $match + $mismatch + $len_mismatch );
	if ( $dbg ) { print "1. match=$match  mismatch=$mismatch  len_mismatch=$len_mismatch  score=$score\n" }
#	if ( uc substr( $peptide, 0, 1 ) eq "M" ) {
#		$len_mismatch -= 5;
#		if ( $len_mismatch < 0 ) { $len_mismatch = 0 }
#		$match += 5;
#		$score = 100.0 * $match / ( $match + $mismatch + $len_mismatch );
#		if ( $dbg ) { print "2. match=$match  mismatch=$mismatch  len_mismatch=$len_mismatch  score=$score\n" }
#	}
	return $score;
}

# return unique elements of an array
sub unique {
	my @list = @_;

	my @uniq;
	my %seen;
	for my $item ( @list ) {
		push(@uniq, $item) unless $seen{$item}++;
	}

	return @uniq;
}

# compare two proteins by aligning them and scoring the alignment
sub quick_seqscore {
	my ( $seq, $ref_id, $refsq ) = @_;

	my $seqlen = length( $seq );

	my $refseq;
	if ( defined $refsq ) {
		$refseq = $refsq;
	}
	else {
		my $ref = get_reference_seq( $ref_id );
		$refseq = $$ref{sequence};
	}
	my $reflen = length( $refseq );
	if ( $refseq =~ /\*$/ ) { $reflen-- }

	my $ratio = 100.0 * $seqlen / $reflen;
	if ( $ratio < 0.75 * get_parameter( "min_generef_pctcoverage" ) ) { return 0 }

	$ratio = 100.0 * $reflen / $seqlen;
	if ( $ratio < 0.75 * get_parameter( "min_generef_pctsimilarity" ) ) { return 0 }

	my $alignment = align_pep_to_reference( "seq", $seq, $ref_id, $refseq );
	my ( $alignlen, $numcov, $numid, $numsim, $num25, $numt5, $numt3 ) = get_alignment_counts( "seq", $ref_id, $alignment );

	my $pctcov = int( 100.0 * $numcov / $reflen + 0.5 );
	my $nummat = ( $numid + $numsim ) / 2.0;
	my $pctsim = int( 100.0 * $numsim  / $alignlen + 0.5 );
	my $numextra = $alignlen - $numcov;
#	my $numins = ( $alignlen - $numcov );
#	my $pctins = int( 100.0 * $numins / $alignlen + 0.5 );
#	my $coverage = $numcov / $reflen;

	if ( $pctsim < 0.75 * get_parameter( "min_generef_pctsimilarity" ) ) { return 0 }

	my $score = 100.0
		* ( $numid/$alignlen + $numsim/$alignlen + $num25/minval(25.0,$alignlen)/4.0  - $numextra/$alignlen/10.0 ) / 2.25
		* $numcov / $reflen;
#	my $score = sqrt( $nummat * $numcov ) - ( $numins );
#print "\t cov: $pctcov\t mat: $pctmat\t ins: $pctins\t score: $score\n";

	return maxval( $score, 0 );
}

# convert a failed gene into a "mutation"
sub gene2mutation {
	my ( $genome, $rawgene, $gene_warning ) = @_;
	my $detailed_warnings = get_parameter( "detailed_exceptions" );

	my $mutation;
	$$mutation{is_mutation} = 1;
	$$mutation{error} = 0;
	if ( defined $gene_warning && $detailed_warnings ) {
		$$mutation{warning} = $gene_warning;
	}
	else {
		$$mutation{warning} = "Possible Sequence Mutations";
	}

	$$mutation{gene_id} = "$$rawgene{query_id}.$$rawgene{gene_num}";

	$$mutation{dna_id}      = $$rawgene{query_id};
	$$mutation{orientation} = $$rawgene{orientation};

	if ( $$rawgene{orientation} == -1 ) {
		$$mutation{start_codon} = $$rawgene{query_right};
		$$mutation{stop_site}   = $$rawgene{query_left};
	}
	else {
		$$mutation{start_codon} = $$rawgene{query_left};
		$$mutation{stop_site}   = $$rawgene{query_right};
	}
	$$mutation{start_truncated} = 0;
	$$mutation{stop_truncated} = 0;

	$$mutation{cdna} = subsequence( $$genome{sequence}, $$mutation{start_codon}, $$mutation{stop_site} );

	my $exon;
	$$exon{dna_start} = $$mutation{start_codon};
	$$exon{dna_end} = $$mutation{stop_site};
	$$exon{cdna_start} = 1;
	$$exon{cdna_end} = $$rawgene{query_right} - $$rawgene{query_left} + 1;
	$$exon{subject_start} = $$rawgene{subject_left};
	$$exon{subject_end} = $$rawgene{subject_right};
	my @exons = ( $exon );
	$$mutation{exons} = \@exons;

	$$mutation{protein} = DNA2AA( $$mutation{cdna} );
	$$mutation{protein_length} = length( $$mutation{protein} );

	$$mutation{ref_id}    = $$rawgene{subject_id};
	$$mutation{ref_start} = $$rawgene{subject_left};
	$$mutation{ref_end}   = $$rawgene{subject_right};
	$$mutation{ref_length}   = $$rawgene{subject_length};
	$$mutation{pct_refcoverage} = $$rawgene{pct_scoverage};
	$$mutation{pct_refidentity} = $$rawgene{pct_identity};
	$$mutation{pct_refsimilarity} = $$rawgene{pct_similarity};

	if ( index( $$mutation{pct_refcoverage}, "." ) < 0 ) {
		$$mutation{pct_refcoverage} .= ".0";
	}
	if ( index( $$mutation{pct_refidentity}, "." ) < 0 ) {
		$$mutation{pct_refidentity} .= ".0";
	}
	if ( index( $$mutation{pct_refsimilarity}, "." ) < 0 ) {
		$$mutation{pct_refsimilarity} .= ".0";
	}

	$$mutation{gene_name} = get_reference_name( $$mutation{ref_id} );
	$$mutation{gene_definition} = get_reference_definition( $$mutation{ref_id} );

	$$mutation{ref_alignment} = bl2seq_alignment( $genome, $mutation, get_parameter( "pep_bl2seqopts") );
	$$mutation{splicing_quality} = 0;
	$$mutation{match_quality} = ( $$mutation{pct_refidentity} + $$mutation{pct_refsimilarity} ) / 200.0
		* $$mutation{pct_refcoverage} / 100.0;
	$$mutation{gene_quality} = $$mutation{match_quality};
	if ( $$rawgene{has_frameshift} ) {
		my $fs_message = "probable frameshift in genome";
		if ( $$mutation{warning} && index( $$mutation{warning}, "%" ) < 0  ) {
			$$mutation{warning} .= ", $fs_message";
		}
		else {
			$$mutation{warning} = $fs_message;
		}
	}

	return $mutation;
}

#
sub calculate_splicing_quality {
	my ( $genome, $gene ) = @_;

	if ( $$gene{is_mutation} ) { return 0 }

	my @exons = @{$$gene{exons}};
	if ( @exons < 2 ) { return 0 }

	my %donorlist = %{get_splice_donors()};
	my $quality = 0;
	my $donor = subsequence( $$genome{sequence}, $exons[0]{dna_end} + $$gene{orientation}, $exons[0]{dna_end} + 2 * $$gene{orientation} );
	for my $i ( 2..@exons ) {
		my $acceptor = subsequence( $$genome{sequence}, $exons[$i-1]{dna_start} - 2 * $$gene{orientation}, $exons[$i-1]{dna_start} - $$gene{orientation} );
		if ( ! defined $donorlist{$donor}{$acceptor} ) {
#print "invalid splice: $donor/$acceptor\n";
			$quality -= 25;
		}
		else {
			$quality -= $donorlist{$donor}{$acceptor};
#print " $donor/$acceptor=$donorlist{$donor}{$acceptor}\n";
		}
		if ( $i < @exons ) {
			$donor = subsequence( $$genome{sequence}, $exons[$i-1]{dna_end} + $$gene{orientation}, $exons[$i-1]{dna_end} + 2 * $$gene{orientation} );
		}
	}
#print "quality=$quality\n";
	return $quality;
}

# calculate quality of alignment between reference and gene
sub score_pep_vs_reference {
	my ( $pep_id, $pep, $ref_id, $ref_pep, $stats ) = @_;

	my $ref_length = length( $ref_pep );
	if ( substr( $ref_pep, $ref_length - 1, 1 ) eq "*" ) {
		$ref_length--;
	}

	my $alignment = align_pep_to_reference( $pep_id, $pep, $ref_id, $ref_pep );
#print "Q $pep_id vs S $ref_id alignment=\n$alignment\n";
	my ( $alignlen, $num_covered, $num_identical, $num_similar, $num_sim25, $num_trunc5, $num_trunc3 ) =
		get_alignment_counts( $pep_id, $ref_id, $alignment );

	( $$stats{alignlen},$$stats{num_covered}, $$stats{num_identical}, $$stats{num_similar},
			$$stats{num_sim25}, $$stats{num_trunc5}, $$stats{num_trunc3} ) =
		( $alignlen, $num_covered, $num_identical, $num_similar, $num_sim25, $num_trunc5, $num_trunc3 );

	if ( $alignlen == 0 ) {
		$$stats{pct_refcoverage} = 0;
		$$stats{pct_refsimilarity} = 0;
		$$stats{pct_refsimilarity25} = 0;
		$$stats{pct_refidentity} = 0;
		$$stats{pct_reftrunc5} = 0;
		$$stats{pct_reftrunc3} = 0;
		$$stats{match_quality} = 0;
	}
	else {
		$$stats{pct_refcoverage} = 100.0 * $num_covered / $ref_length;
		$$stats{pct_refsimilarity} = 100.0 * $num_similar / $alignlen;
		$$stats{pct_refsimilarity25} = 100.0 * $num_sim25 / minval( 25.0, $alignlen );
		$$stats{pct_refidentity} = 100.0 * $num_identical / $alignlen;
		$$stats{pct_reftrunc5} = 100.0 * $num_trunc5 / $ref_length;
		$$stats{pct_reftrunc3} = 100.0 * $num_trunc3 / $ref_length;
		$$stats{match_quality} = 100.0 *
			( $num_identical / $alignlen
				+ $num_similar / $alignlen
				+ $num_sim25 / minval( 25.0, $alignlen ) / 4.0 ) / 2.25
			* $num_covered / $ref_length;
	}

	$$stats{pep_alignment} = $alignment;
}

sub matpep_alignment {
	my ( $gene ) = @_;

	if ( ! defined $$gene{mature_peps} ) { return }

	my $refs = 0;
	my $mpseq = "";
	for my $pep ( @{$$gene{mature_peps}} ) {
		my $seq = $$pep{pep_sequence};
		if ( $$pep{ref_id} eq "?" ) {
			$seq =~ s/./X/g;
		}
		else {
			my $ref = $refmat{$$pep{ref_id}};
			$seq = uc $$ref{sequence};
			$refs++;
		}
		if ( length( $mpseq ) ) { $mpseq .= "xxxxx" }
		$mpseq .= $seq;
	}
	if ( $refs < 2 ) { return }
	my $alignment = align_pep_to_reference( "MATPEPS", $mpseq, $$gene{gene_id}, $$gene{protein} );
	$$gene{mp_alignment} = $alignment;
}

sub round_refstats {
	my ( $stats ) = @_;
	for my $statname ( get_refstat_attributes() ) {
		if ( $statname eq "alignment" ) { next }
		$$stats{$statname} = int( 10.0 * $$stats{$statname} + 0.5 ) / 10.0;
		if ( index( $$stats{$statname}, "." ) < 0 ) { $$stats{$statname} .= ".0" }
	}
}

sub get_alignment_counts {
	my ( $gene_id, $ref_id, $alignment ) = @_;
#print "gene $gene_id ref $ref_id align $alignment\n";
# parse the genome:reference alignment
	my $dbg = 0;
	#$dbg = 1;

	my @tmp = split /\n/, $alignment;
	shift @tmp;
	shift @tmp;
	shift @tmp;

	my $i        = 0;
	my $genestr  = "";
	my $refstr   = "";
	my $compstr  = "";
	my $prevstr = "";
	my $margin;
	my $field = 0;
	my $genecheck = $gene_id;
	if ( length( $genecheck ) > 30 ) { $genecheck = substr( $genecheck, 0, 30 ) }
	for my $line ( @tmp ) {
		my $seqid = $line;
		$seqid =~ s/ .*$//;
		if ( $seqid eq $genecheck ) {
			$field = 1;
			my $alignstr = $line;
			$alignstr =~ s/^[^ ]*  *//;
			$genestr .= $alignstr;
			$margin = index( $line, $alignstr );
			$prevstr = $alignstr;
		}
		elsif ( $seqid eq $ref_id ) {
			$field = 2;
			my $alignstr = $line;
			$alignstr =~ s/^[^ ]*  *//;
			$refstr .= $alignstr;
			$margin = index( $line, $alignstr );
			$prevstr = $alignstr;
		}
		else {
			$field++;
			if ( $field == 3 ) {
				if ( $line !~ /^ *$/ ) {
					my $alignstr = substr( $line, $margin );
					$compstr .= $alignstr;
				}
				else {
					my $alignstr = $prevstr;
					$alignstr =~ s/[^ ]/ /g;
					$compstr .= $alignstr;
				}
			}
			else {
				$field = 4;
				next;
			}
		}
	}
	if ( $dbg ) { print "\nalignment\ngene=$genestr\ncomp=$compstr\n ref=$refstr\n" }
	my $tmp = $refstr;
	$tmp =~ s/-//g;
	my $reflen = length( $tmp );

	while ( length( $compstr ) < length( $refstr ) ) {
		$compstr .= "          ";
	}
	if ( $dbg ) { print "reflen=$reflen\n" }

# check for 5' and 3' truncation
	my $num_trunc5 = 0;
	if ( $genestr =~ /^(\-+)/ ) {
		$num_trunc5 = length($1);
		$genestr = substr( $genestr, $num_trunc5 );
		$compstr = substr( $compstr, $num_trunc5 );
		$refstr = substr( $refstr, $num_trunc5 );
	}
	my $num_trunc3 = 0;
	if ( $genestr =~ /(\-+)$/ ) {
		$num_trunc3 = length($1);
		$genestr = substr( $genestr, 0, length($genestr) - $num_trunc3 );
		$compstr = substr( $compstr, 0, length($compstr) - $num_trunc3 );
		$refstr = substr( $refstr, 0, length($refstr) - $num_trunc3 );
	}
#	my $alignlen = length( $refstr );
	if ( $dbg ) { print "\nafter truncation\ngene=$genestr\ncomp=$compstr\n ref=$refstr\n" }

# remove internal gaps on gene
	my $pos = 0;
	my $tmpgene = $genestr;
	$genestr = "";
	my $tmpcomp = $compstr;
	$compstr = "";
	my $tmpref = $refstr;
	$refstr = "";
	for my $frag ( split /([xX-]{4,})/, $tmpgene ) {
		if ( $frag !~ /[xX-]{4,}/ ) {
			$genestr .= substr( $tmpgene, $pos, length( $frag ) );
			$compstr .= substr( $tmpcomp, $pos, length( $frag ) );
			$refstr .= substr( $tmpref, $pos, length( $frag ) );
		}
		$pos += length( $frag )
	}
	if ( $dbg ) { print "\nafter gene gaps\ngene=$genestr\ncomp=$compstr\n ref=$refstr\n" }
	my $alignlen = length( $refstr );

# calculate coverage and matches
	$refstr =~ s/-//g;
	my $num_covered = length( $refstr );

	my $comp25 = substr( $compstr, 0, 25 );
	$comp25 =~ s/ //g;
	my $num_similar25 = length( $comp25 );
	$comp25 =~ s/\.//g;
	my $lowsim = $num_similar25 - length( $comp25 );
	$num_similar25 -= $lowsim / 4.0;

	$compstr =~ s/ //g;
	my $num_similar = length($compstr);
	$compstr =~ s/\.//g;
	$lowsim = $num_similar - length( $compstr );
	$num_similar -= $lowsim / 4.0;

	$compstr =~ s/[^*]//g;
	my $num_identical = length($compstr);

	if ( $dbg ) {
		print "alignlen=$alignlen  num_trunc5=$num_trunc5  num_trunc3=$num_trunc3\n";
		print "num_covered=$num_covered  num_similar=$num_similar  num_identical=$num_identical\n";
	}

	return ( $alignlen, $num_covered, $num_identical, $num_similar, $num_similar25,
		$num_trunc5, $num_trunc3 );
}

# convert peptide sequence to dna regular expression
# ordered: 0=allow AAs in any order, 1=require AAs match in order
sub pep_to_dnaexp {
	my ( $pep, $ordered ) = @_;
	if ( ! defined $ordered ) { $ordered = 1 }

	my %codons = %{ get_codons() };
	my $dna = "";

	for my $i ( 0 .. length( $pep ) - 1 ) {
		my $aa = uc substr( $pep, $i, 1 );
		my @tmp;
		my @exps;
		for my $codon ( keys %codons ) {
			if ( uc $codons{$codon} eq $aa ) {
				push @exps, $codon;
			}
		}
		if ( ! $ordered ) {
			$dna .= join( "|", @exps );
		}
		elsif ( @exps > 1 ) {
			$dna .= "(" . join( "|", @exps ) . ")";
		}
		elsif ( index( $exps[0], "|" ) >= 0 ) {
			$dna .= "($exps[0])";
		}
		else {
			$dna .= $exps[0];
		}
	}

	return $dna;
}

sub validate_reference_alignment {
	my ( $genome, $gene ) = @_;

	my $validation = before_validate_gene( $genome, $gene );
	if ( defined $validation ) { return $validation }

	my $min_coverage = get_parameter("min_generef_pctcoverage");
	if ( ! check_coverage( $min_coverage, $gene, $genome, 0 ) ) {
		$$gene{is_mutation} = 1;
		my $msg = "cds covers $$gene{pct_refcoverage}% of reference, at least $min_coverage% required";
		if ( defined $$gene{warning} && length( $$gene{warning} ) ) {
			$$gene{warning} .= ", $msg";
		}
		else {
			$$gene{warning} = $msg;
		}
		$validation = after_validate_gene( $genome, $gene );
		if ( defined $validation ) { return $validation }
		return 0;
	}

	my $min_similarity = get_parameter("min_generef_pctsimilarity");
	if ( $$gene{pct_refsimilarity} < $min_similarity ) {
		$$gene{is_mutation} = 1;
		my $msg = "cds is $$gene{pct_refsimilarity}% similar to reference, at least $min_similarity% required";
		if ( defined $$gene{warning} && length( $$gene{warning} ) ) {
			$$gene{warning} .= ", $msg";
		}
		else {
			$$gene{warning} = $msg;
		}
		$validation = after_validate_gene( $genome, $gene );
		if ( defined $validation ) { return $validation }
		return 0;
	}

	my $min_similarity25 = get_parameter("min_generef_pctsimilarity25");
	if ( $$gene{pct_refsimilarity25} + $$gene{pct_reftrunc5} < $min_similarity25 ) {
		$$gene{is_mutation} = 1;
		my $msg = "cds is $$gene{pct_refsimilarity25}% similar to reference over first 25 aa, at least $min_similarity25% required";
		if ( defined $$gene{warning} && length( $$gene{warning} ) ) {
			$$gene{warning} .= ", $msg";
		}
		else {
			$$gene{warning} = $msg;
		}
		$validation = after_validate_gene( $genome, $gene );
		if ( defined $validation ) { return $validation }
		return 0;
	}

	my $min_identity = get_parameter("min_generef_pctidentity");
	if ( $$gene{pct_refidentity} < $min_identity ) {
		$$gene{is_mutation} = 1;
		my $msg = "cds is $$gene{pct_refidentity}% identical to reference, at least $min_identity% required";
		if ( defined $$gene{warning} && length( $$gene{warning} ) ) {
			$$gene{warning} .= ", $msg";
		}
		else {
			$$gene{warning} = $msg;
		}
		$validation = after_validate_gene( $genome, $gene );
		if ( defined $validation ) { return $validation }
		return 0;
	}

	if ( ! $$gene{start_truncated} ) {
		my $max_trunc5 = get_parameter("max_generef_pcttrunc5");
		if ( $$gene{pct_reftrunc5} > $max_trunc5 ) {
			$$gene{is_mutation} = 1;
			my $msg = "cds 5' truncation of $$gene{pct_reftrunc5}% exceeds maximum allowed value: $max_trunc5%";
			if ( defined $$gene{warning} && length( $$gene{warning} ) ) {
				$$gene{warning} .= ", $msg";
			}
			else {
				$$gene{warning} = $msg;
			}
			$validation = after_validate_gene( $genome, $gene );
			if ( defined $validation ) { return $validation }
			return 0;
		}
	}
	if ( ! $$gene{stop_truncated} ) {
		my $max_trunc3 = get_parameter("max_generef_pcttrunc3");
		if ( $$gene{pct_reftrunc3} > $max_trunc3 ) {
			$$gene{is_mutation} = 1;
			my $msg = "cds 3' truncation of $$gene{pct_reftrunc3}% exceeds maximum allowed value: $max_trunc3%";
			if ( defined $$gene{warning} && length( $$gene{warning} ) ) {
				$$gene{warning} .= ", $msg";
			}
			else {
				$$gene{warning} = $msg;
			}
			$validation = after_validate_gene( $genome, $gene );
			if ( defined $validation ) { return $validation }
			return 0;
		}
	}

	$validation = after_validate_gene( $genome, undef, $gene );
	if ( defined $validation ) { return $validation }
	return 1;
}

sub mapto_polyprotein {
	my ( $gene, $polyprotein, $matpepdb, $cleavagesites ) = @_;

	if ( $$gene{is_mutation} ) { return }

	if ( ! defined $matpepdb ) {
		$matpepdb = get_parameter( "mature_pep_refdb" );
	}

	my $matpep_blastopts = get_parameter( "mature_pep_blastopts" );
	my $matpep_mincoverage = get_parameter( "mature_pep_mincoverage" );
	my $matpep_minsimilarity = get_parameter( "mature_pep_minsimilarity" );
	my $allow_fuzzy_matpeps = get_parameter( "allow_fuzzy_matpeps" );

	if ( !defined $cleavagesites ) {
		$cleavagesites = get_cleavage_sites();
	}

	if ( ! defined $polyprotein ) { $polyprotein = $$gene{protein} };
	my $polylen = length( $polyprotein );
	if ( substr( $polyprotein, $polylen-1, 1 ) eq "*" ) { $polylen-- }

	my @rawmatpeps;
	{
		my @tmp = sort { $$a{query_left} <=> $$b{query_left} }
			blast_sequence( "$$gene{gene_id}", $polyprotein, $matpepdb,
				$matpep_blastopts, 0 );

		for my $hsp ( @tmp ) {
			if ( ! check_coverage( $matpep_mincoverage, $hsp, $gene, 0 ) ) { next }


			# clean up reference definition
			if ( $$hsp{subject_definition} =~ /^([^ ]+) .*product="([^"]+)"/ ) {
				$$hsp{subject_definition} = "$1 $2";
			}

			# extrapolate edges to full coverage of subject sequence
			if ( $$hsp{pct_similarity} < $matpep_minsimilarity ) { next }
			my $proj_left = $$hsp{query_left} - ( $$hsp{subject_left} - 1 );
			my $proj_right = $$hsp{query_right} + ( $$hsp{subject_length} - $$hsp{subject_right} );
			if ( $proj_left < 1 ) { $proj_left = 1 }
			if ( $proj_right > $polylen ) { $proj_right = $polylen }
			if ( $proj_left < $$hsp{query_left} ) { $$hsp{vigor_matchwgt} += $proj_left - $$hsp{query_left} }
			if ( $proj_right > $$hsp{query_right} ) { $$hsp{vigor_matchwgt} += $$hsp{query_right} - $proj_right }

			$$hsp{proj_left} = $proj_left;
			$$hsp{proj_right} = $proj_right;
			push @rawmatpeps, $hsp;

			if ( $$hsp{proj_left} != $$hsp{query_left} ) {
				my %tmp = %$hsp;
				( $tmp{proj_left}, $tmp{query_left} )   = ( $tmp{query_left}, $tmp{proj_left} );
				push @rawmatpeps, \%tmp;

				if ( $tmp{proj_right} != $tmp{query_right} ) {
					my %tmp2 = %tmp;
					( $tmp2{proj_right}, $tmp2{query_right} ) = ( $tmp2{query_right}, $tmp2{proj_right} );
					push @rawmatpeps, \%tmp2;
				}
			}
			elsif ( $$hsp{proj_right} != $$hsp{query_right} ) {
				my %tmp = %$hsp;
				( $tmp{proj_right}, $tmp{query_right} ) = ( $tmp{query_right}, $tmp{proj_right} );
				push @rawmatpeps, \%tmp;
			}
		}
	}
	if ( ! @rawmatpeps ) { return undef }

	# eliminate mature peptides with identical edges
	# keeping peptide with best similarity
	my $i = 0;
	while ( $i < @rawmatpeps ) {

		my $ihit = $rawmatpeps[$i];
		my $j = $i + 1;
		while ( $j < @rawmatpeps ) {
			my $jhit = $rawmatpeps[$j];

			if ( $$ihit{proj_left}==$$jhit{proj_left} && $$ihit{proj_right}==$$jhit{proj_right} ) {
				if ( $$ihit{vigor_pctsimilarity} >= $$jhit{vigor_pctsimilarity} ) {
					splice @rawmatpeps, $j, 1;
					next;
				}
				else {
					splice @rawmatpeps, $i, 1;
					$i--;
					last;
				}
			}
			else {
				$j++;
			}
		}
		$i++;
	}

	# score the edges based on similarity of matches
	my %lefts;
	my %rights;
	my $errbar = 5;
	my $errval = 2.5;
	{
		my $pct = 100;
		my $base = 1;
		my $edge = \%lefts;
		for my $i ( $base-$errbar .. $base+$errbar ) {
			$$edge{$i}{similarity} = $pct - $errval * abs( $i - $base );
		}
	}
	{
		my $pct = 100;
		my $base = $polylen;
		my $edge = \%rights;
		for my $i ( $base-$errbar .. $base+$errbar ) {
			$$edge{$i}{similarity} = $pct - $errval * abs( $i - $base );
		}
	}
	for my $hsp ( sort { $$b{vigor_pctsimilarity} <=> $$a{vigor_pctsimilarity} } @rawmatpeps ) {
		my $pct = $$hsp{vigor_pctsimilarity};
		if ( $pct < 40 ) { $pct = 40 }
		my $base = $$hsp{proj_left}-1;
		my $edge = \%rights;
		for my $i ( $base-$errbar .. $base+$errbar ) {
			$$edge{$i}{similarity} = $pct - $errval * abs( $i - $base );
		}
		$base = $$hsp{proj_right}+1;
		$edge = \%lefts;
		for my $i ( $base-$errbar .. $base+$errbar ) {
			$$edge{$i}{similarity} = $pct - $errval * abs( $i - $base );
		}
	}

	# score edges on presence of likely cleavage sites
	for my $edge ( sort { $a <=> $b } keys %lefts ) {
		$lefts{$edge}{cleavage} = 0;
		if ( defined $cleavagesites ) {
			for my $regexp ( sort { $$cleavagesites{$b}{bonus}
										<=> $$cleavagesites{$b}{bonus} } keys %$cleavagesites ) {
				my $offset = $$cleavagesites{$regexp}{offset};
				my $start = $edge - $offset - 1;
				if ( $start < 0 ) { next }
				my $len = minval( $polylen - $start, 10 );
				if ( $len < 1 ) { next  }
				my $seq = substr( $polyprotein, $start, $len );
				if ( $seq =~ /^$regexp/ ) {
					$lefts{$edge}{cleavage} = $$cleavagesites{$regexp}{bonus};
					last;
				}
			}
		}
	}

	for my $edge ( sort { $a <=> $b } keys %rights ) {
		$rights{$edge}{cleavage} = 0;
		if ( defined $cleavagesites ) {
			for my $regexp ( sort { $$cleavagesites{$b}{bonus}
										<=> $$cleavagesites{$b}{bonus} } keys %$cleavagesites ) {
				my $offset = $$cleavagesites{$regexp}{offset};
				my $start = $edge - $offset;
				if ( $start < 0 ) { next }
				my $len = minval( $polylen - $start, 10 );
				if ( $len < 1 ) { next }
				my $seq = substr( $polyprotein, $start, $len );
				if ( $seq =~ /^$regexp/ ) {
					$rights{$edge}{cleavage} = $$cleavagesites{$regexp}{bonus};
					last;
				}
			}
		}
	}

	# score the mature peptides based on their edges
	@rawmatpeps = sort { $$a{proj_left} <=> $$b{proj_left} } @rawmatpeps;
	for my $hsp ( @rawmatpeps ) {
		$$hsp{edge_score} = 0.5 * $$hsp{vigor_pctsimilarity};
		if ( exists $lefts{$$hsp{proj_left}} ) {
			$$hsp{edge_score} += 0.25 * $lefts{$$hsp{proj_left}}{similarity}
				+ $lefts{$$hsp{proj_left}}{cleavage};
		}
		if ( exists $rights{$$hsp{proj_right}} ) {
			$$hsp{edge_score} += 0.25 * $rights{$$hsp{proj_right}}{similarity}
				+ $rights{$$hsp{proj_right}}{cleavage};
		}
		$$hsp{vigor_matchwgt} = $$hsp{edge_score};
	}

	if ( get_parameter( "verbose" ) ) {
		print "\nSCORED MATURE PEPS (see vwgt)\n";
		print_blasthits( @rawmatpeps );
	}

	# eliminate overlapping mature peptides
	# keep best scoring peptide
	$i = 0;
	my @saved;
	while ( $i < @rawmatpeps ) {
		my $ihit = $rawmatpeps[$i];
		$$ihit{start_fuzzy} = 0;
		$$ihit{end_fuzzy} = 0;

		my $j = $i + 1;
		while ( $j < @rawmatpeps ) {
			my $jhit = $rawmatpeps[$j];

			my $overlap = minval( $$ihit{proj_right}, $$jhit{proj_right} )
				- maxval(  $$ihit{proj_left},  $$jhit{proj_left} );
			my $max_allowed = minval( $$ihit{subject_coverage}, $$jhit{subject_coverage} ) / 3.0;
			if  ( $max_allowed < 5 ) { $max_allowed = 5 }
			if ( $overlap > $max_allowed ) {
				if ( $$ihit{edge_score} + $$ihit{vigor_pctsimilarity} >
						$$jhit{edge_score} + $$jhit{vigor_pctsimilarity}) {
					push @saved, splice @rawmatpeps, $j, 1;
					next;
				}
				elsif ( $$ihit{edge_score} < $$jhit{edge_score} ) {
					push @saved, splice @rawmatpeps, $i, 1;
					$i--;
					last;
				}
				elsif ( $$ihit{vigor_pctsimilarity} >= $$jhit{vigor_pctsimilarity} ) {
					push @saved, splice @rawmatpeps, $j, 1;
					next;
				}
				else {
					push @saved, splice @rawmatpeps, $i, 1;
					$i--;
					last;
				}
			}
			else {
				$j++;
			}
		}
		$i++;
	}

	if ( get_parameter( "verbose" ) ) {
		print "\nSELECTED MATURE PEPS\n";
		print_blasthits( @rawmatpeps );
	}

	# adjust edges to eliminate gaps
	%refmat = loadFasta( $matpepdb );
	{
		my @tmp = sort { $$a{proj_left} <=> $$b{proj_left} } @rawmatpeps;
		@rawmatpeps = ();

		# make sure mature peptides start at first aa
		my $right = 0;
		{
			my $raw = $tmp[0];

			# large gap implies a missing mature peptide
			my $maxgap = maxval( 10, $$raw{subject_length} / 10.0 );
			if ( $$raw{proj_left} > $maxgap ) {

				my $new = best_matpep_match( 1, $right, @saved );
				if ( ! defined $new ) {
					my %tmp = %$raw;
					$new = \%tmp;
					$$new{subject_id} = "?";
					$$new{subject_definition} = "? unknown mature peptide";
				}
				$$new{proj_left} = 1;
				$$new{proj_right} = $$raw{proj_left}-1;
				( $$new{query_left} , $$new{query_right} ) = ( $$new{proj_left}, $$new{proj_right}-2 );
				$$new{start_fuzzy} = 0;
				$$new{end_fuzzy} = 0;
				adjust_edge( $polyprotein, \%refmat, $new, $raw, $cleavagesites, $allow_fuzzy_matpeps );

				push @rawmatpeps, $new;
				$right = $$new{proj_right};
			}
			else {
				$$raw{proj_left} = 1;
			}
		}

		# mature peptides should be adjacent
		# adjust alignment boundaries appropriately
		while ( @tmp > 1 ) {
			my $raw = shift @tmp;
			push @rawmatpeps, $raw;

			if ( ! $allow_fuzzy_matpeps ) {
				$$raw{proj_left} = $right + 1;
			}

			# large gap implies a missing mature peptide
			my $rawnext = $tmp[0];
			my $maxgap = maxval( 10, maxval( $$rawnext{subject_length}, $$raw{subject_length} ) / 10.0 );
			my $gap = $$rawnext{proj_left} - $$raw{proj_right};
			if ( $gap > $maxgap ) {
				my $new = best_matpep_match( $$raw{proj_right}+1, $$rawnext{proj_left}-1, @saved );
				if ( ! defined $new ) {
					my %tmp = %$raw;
					$new = \%tmp;
					$$new{subject_id} = "?";
					$$new{subject_definition} = "? unknown mature peptide";
				}
				$$new{proj_left} = $$raw{proj_right}+1;
				$$new{proj_right} = $$rawnext{proj_left}-1;
				( $$new{query_left} , $$new{query_right} ) = ( $$new{proj_left}+2, $$new{proj_right}-2 );
				$$new{start_fuzzy} = 0;
				$$new{end_fuzzy} = 0;
				adjust_edge( $polyprotein, \%refmat, $raw, $new, $cleavagesites, $allow_fuzzy_matpeps );
				adjust_edge( $polyprotein, \%refmat, $new, $rawnext, $cleavagesites, $allow_fuzzy_matpeps );

				push @rawmatpeps, $new;
				$raw = $new;
			}
			adjust_edge( $polyprotein, \%refmat, $raw, $rawnext, $cleavagesites, $allow_fuzzy_matpeps );
		}

		# make sure mature peptides end at last aa
		{
			my $raw = $tmp[0];
			push @rawmatpeps, $raw;

			if ( ! $allow_fuzzy_matpeps ) {
				$$raw{proj_left} = $right + 1;
			}

			# large gap implies a missing mature peptide
			my $maxgap = maxval( 10, $$raw{subject_length} / 10.0 );
			if ( $polylen+1 - $$raw{proj_right} > $maxgap ) {
				my $new = best_matpep_match( $right+1, $polylen, @saved );
				if ( ! defined $new ) {
					my %tmp = %$raw;
					$new = \%tmp;
					$$new{subject_id} = "?";
					$$new{subject_definition} = "? unknown mature peptide";
				}
				$$new{proj_left} = $$raw{proj_right}+1;
				$$new{proj_right} = $polylen;
				( $$new{query_left} , $$new{query_right} ) = ( $$new{proj_left}+2, $$new{proj_right} );
				$$new{start_fuzzy} = 0;
				$$new{end_fuzzy} = 0;
				adjust_edge( $polyprotein, \%refmat, $raw, $new, $cleavagesites, $allow_fuzzy_matpeps );

				push @rawmatpeps, $new;
			}
			else {
				$$raw{proj_right} = $polylen;
			}
		}
	}

	if ( get_parameter( "verbose" ) ) {
		print "\nPOLISHED MATURE PEPS\n";
		print_blasthits( @rawmatpeps );
	}

	# return the final list of mature peptides
	my $pepno = 0;
	my @matpeps;
	for my $raw ( @rawmatpeps ) {
		$pepno++;

		my $pep;
		$$pep{pep_id}         = "$$gene{gene_id}.$pepno";
		$$pep{pep_definition} = $$raw{subject_definition};
		$$pep{pep_definition} =~ s/^[^ ]* *//;
		$$pep{pep_start} = $$raw{proj_left};
		$$pep{pep_end} = $$raw{proj_right};
		$$pep{pep_length} = $$pep{pep_end}-$$pep{pep_start}+1;
		$$pep{pep_sequence} = substr( $polyprotein, $$pep{pep_start}-1, $$pep{pep_length} );
		$$pep{start_fuzzy} = $$raw{start_fuzzy};
		$$pep{end_fuzzy} = $$raw{end_fuzzy};
		$$pep{ref_id} = $$raw{subject_id};
		if ( $$raw{subject_id} ne "?" ) {
			my $refpep = $refmat{ $$pep{ref_id} }{sequence};
			score_pep_vs_reference( $$pep{pep_id}, $$pep{pep_sequence}, $$pep{ref_id}, $refpep, $pep );
			delete( $$pep{pep_alignment} );
			round_refstats( $pep );
		}

		push @matpeps, $pep;
	}
	if ( ! $matpeps[0]{start_fuzzy} ) {
		$matpeps[0]{start_fuzzy} = $$gene{start_truncated};
	}
	if ( ! $matpeps[@matpeps-1]{end_fuzzy} ) {
		$matpeps[@matpeps-1]{end_fuzzy} = $$gene{stop_truncated};
	}

	return \@matpeps;
}

sub best_matpep_match {
	my ( $left, $right, @candidates ) = @_;
	my $match;
	my $matchscore;
#print "find best fit for $left-$right\n";
	for my $candidate ( sort { compare_qry_positions( $a, $b ) } @candidates ) {
		if ( $$candidate{query_right} <  $left ) { next }
		if ( $$candidate{query_left} >  $right ) { last }
#print_blasthits( $candidate );
		my $overlap =
			( minval( $right, $$candidate{query_right} )
				- maxval(  $left,  $$candidate{query_left} ) );
		my $region_coverage = int( 100.0 * $overlap / ( $right-$left+1 ) + 0.5 );
		my $matpep_coverage = int( 100.0 * $overlap / $$candidate{subject_coverage} + 0.5 );
#print "  overlap=$overlap  region coverage=$region_coverage  matpep coverage=$matpep_coverage\n";
		if ( $region_coverage < 66.0 ) { next }
		if ( $matpep_coverage < 67.0 ) { next }
		my $match_score = $overlap * $$candidate{pct_similarity} / 100.0;
		my $mismatch_score =
			( $right - $left + 1
				+ $$candidate{query_right} - $$candidate{query_left} + 1 )
			* 0.25;
		my $score = $match_score - $mismatch_score;
#print " MATCH score=$score";
		if ( ! defined $matchscore || $score > $matchscore ) {
			$matchscore = $score;
			$match = $candidate;
#print " ***";
		}
#print "\n";
	}

	return $match;
}

sub adjust_edge {
	my ( $polyprotein, $refmat, $mpleft, $mpright, $cleavagesites, $allow_fuzzy_edges ) = @_;

	my $dbg = 0;

	my $method = 0;
	my $polylen = length( $polyprotein );
	if ( $polyprotein =~ /\*$/ ) { $polylen-- }

	my $lenl = 0;
	my $refl;
	my %profl;
	if ( $$mpleft{subject_id} ne "?" ) {
		$refl = $$refmat{$$mpleft{subject_id}}{sequence};
		if ( $method == 0 ) { %profl = profile_peptide( $refl ) }
		$lenl = length( $refl );
	}
	my $lenr;
	my $refr;
	my %profr;
	if ( $$mpright{subject_id} ne "?" ) {
		$refr = $$refmat{$$mpright{subject_id}}{sequence};
		if ( $method == 0 ) { %profr = profile_peptide( $refr ) }
		$lenr = length( $refr );
	}
	if ( ! $lenl && ! $lenr ) {
		if ( $dbg ) { print "no reference\n" }
		return;
	}

	my $bestedge;
	my $bestscore = 0;
	my $bestsite = "";
	my $lowscore = 999999999;
	my $start = minval( $$mpleft{proj_right}, $$mpright{proj_left}-1, $$mpleft{query_right}, $$mpright{query_left}-1 ) - 1;
	my $end = maxval( $$mpleft{proj_right}, $$mpright{proj_left}-1, $$mpleft{query_right}, $$mpright{query_left}-1 ) + 1;
	if ( $$mpleft{subject_id} eq "?" ) { $start -= 3 }
	if ( $$mpright{subject_id} eq "?" ) { $end += 3 }
	if ( $start < 1 ) { $start = 1 }
	if ( $end > $polylen ) { $end = $polylen }

#print "\nL: $$mpleft{subject_id}  R: $$mpright{subject_id}  segment: $$mpleft{proj_left}-$$mpright{proj_right}  edge=$$mpleft{proj_right}|$$mpright{proj_left} $start-$end\n";
	if ( $dbg ) { print "L: $$mpleft{subject_id}  R: $$mpright{subject_id}  segment: $$mpleft{proj_left}-$$mpright{proj_right}  edge=$$mpleft{proj_right}|$$mpright{proj_left} $start-$end\n" }
	for my $edge ( $start..$end ) {

		my $pepscore = 0;
		if ( $lenl ) {
			my $mp = substr( $polyprotein, $$mpleft{proj_left}-1, $edge-$$mpleft{proj_left}+1 );
			if ( $method == 0 ) {
				$pepscore += $lenl * score_profile( $mp, %profl ) / 100.0;
			}
			else {
				$pepscore += $lenl * quick_seqscore( $mp, $$mpleft{subject_id}, $refl ) / 100.0;
			}
		}
		if ( $lenr ) {
			my $mp = substr( $polyprotein, $edge, $$mpright{proj_right}-$edge );
			if ( $method == 0 ) {
				$pepscore += $lenr * score_profile( $mp, %profr ) / 100.0;
			}
			else {
				$pepscore += $lenr * quick_seqscore( $mp, $$mpright{subject_id}, $refr ) / 100.0;
			}
		}

		my $sitescore = 0;
		my $sitekey = "";
		if ( defined $cleavagesites ) {
			for my $site ( sort { $$cleavagesites{$b}{bonus} <=> $$cleavagesites{$a}{bonus} } keys %$cleavagesites ) {
				my $offset = $$cleavagesites{$site}{offset};
				my $pos = $edge - $offset;
				if ( $pos < 0 ) { next }
				if ( substr( $polyprotein, $pos, length( $site ) ) eq $site ) {
					$sitescore = $$cleavagesites{$site}{bonus};
					$sitekey = $site;
					last;
				}
			}
		}

		my $score = $pepscore + $sitescore;
		if ( $score > $bestscore ) {
			$bestscore = $score;
			$bestsite = $sitekey;
			$bestedge = $edge;
		}
		elsif ( $score < $lowscore ) {
			$lowscore = $score;
		}
		if ( $dbg ) { print "  edge=$edge  pepscore=$pepscore  sitescore=$sitescore  range=$lowscore-$bestscore\n" }
#print "  edge=$edge  pepscore=$pepscore  site=$sitekey  sitescore=$sitescore  range=$lowscore-$bestscore\n";
	}
#print "  bestedge=$bestedge  bestcore=$bestscore  bestsite=$bestsite\n";
	if ( $dbg ) { print "  bestedge=$bestedge  bestsite=$bestsite  scorerange=$lowscore-$bestscore\n" }
	if ( $bestsite ne "" ) {
		$$mpleft{proj_right} = $bestedge;
		$$mpright{proj_left} = $bestedge+1;
		if ( $dbg ) { print "return 1. $$mpleft{proj_right}|$$mpright{proj_left}\n" }
	}
	else {
		if ( $dbg ) { print "$$mpleft{subject_id} P $$mpleft{proj_right} Q $$mpleft{query_right}\n" }
		if ( $dbg ) { print "$$mpright{subject_id} P $$mpright{proj_left} Q $$mpright{query_left}\n" }
		if ( $$mpleft{subject_id} ne "?" && $$mpright{subject_id} ne "?" ) {
			if ( $$mpleft{proj_right} == $$mpright{proj_left}-1 ) {
				if ( $dbg ) { print "return 2. $$mpleft{proj_right}|$$mpright{proj_left}\n" }
				return;
			}
			if ( $$mpleft{proj_right} == $$mpright{query_left}-1 ) {
				$$mpright{proj_left} = $$mpright{query_left};
				if ( $dbg ) { print "return 3. $$mpleft{proj_right}|$$mpright{proj_left}\n" }
				return;
			}
			if ( $$mpleft{query_right} == $$mpright{proj_left}-1 ) {
				$$mpleft{proj_right} = $$mpleft{query_right};
				if ( $dbg ) { print "return 4. $$mpleft{proj_right}|$$mpright{proj_left}\n" }
				return;
			}
			if ( $$mpleft{query_right} eq $$mpright{query_left}-1 ) {
				$$mpright{proj_left} = $$mpright{query_left};
				$$mpleft{proj_right} = $$mpleft{query_right};
				if ( $dbg ) { print "return 5. $$mpleft{proj_right}|$$mpright{proj_left}\n" }
				return;
			}
		}
		if ( $$mpleft{subject_id} eq "?" ) {
			$$mpleft{proj_right} = $$mpright{proj_left}-1;
			if ( $dbg ) { print "return 6. $$mpleft{proj_right}|$$mpright{proj_left}\n" }
		}
		elsif ( $$mpright{subject_id} eq "?" ) {
			$$mpright{proj_left} = $$mpleft{proj_right}+1;
			if ( $dbg ) { print "return 7. $$mpleft{proj_right}|$$mpright{proj_left}\n" }
		}
		else {
			$$mpleft{proj_right} = $bestedge;
			$$mpright{proj_left} = $bestedge+1;
			if ( $dbg ) { print "return 8. $$mpleft{proj_right}|$$mpright{proj_left}\n" }
		}
		if ( $allow_fuzzy_edges ) {
			$$mpleft{end_fuzzy} =  1;
			$$mpright{start_fuzzy} = 1;
		}
	}
	return;
}


sub split_polyprotein {
	my ( $gene_id, $poly_protein, $firstpepno, $split_sites, $lastdef, $poly_skip ) = @_;

	if ( ! defined $poly_skip ) { $poly_skip = 0 }
	my $pepno    = $firstpepno;
	my $poly_len = length($poly_protein);
	my @cleavages = sort { compare_cleavage_sites( $a, $b ) } @$split_sites;

	my @matpeps;
	my $poly_pos = 0;
	for my $i ( 0 .. @cleavages - 1 ) {
		my ( $bases, $offset, $pepdef ) = split / *\| */, $$split_sites[$i];
		$offset -= $poly_skip;
		if ( $offset >=0 ) {
			my $pep_split =
			  index( $poly_protein, $bases, $offset ) + length($bases);
			if ( $pep_split < length($bases) ) { $pep_split = $poly_len }

			my $pep;
			$$pep{pep_id}         = "$gene_id.$pepno";
			$$pep{pep_definition} = $pepdef;
			$$pep{pep_sequence} =
			  substr( $poly_protein, $poly_pos, $pep_split - $poly_pos );
			$$pep{pep_length} = length( $$pep{pep_sequence} );
			$$pep{pep_start}  = $poly_pos + 1;
			$$pep{pep_end}    = $pep_split;
			push @matpeps, $pep;

			$poly_pos = $pep_split;
			if ( $poly_pos >= $poly_len ) { last }
		}

		$pepno++;
	}

	# get last peptide
	# (remainder after last split)
	if ( defined $lastdef && $poly_pos < $poly_len ) {
		my $pep;
		$$pep{pep_id}         = "$gene_id.$pepno";
		$$pep{pep_definition} = $lastdef;
		$$pep{pep_start}      = $poly_pos + 1;
		$$pep{pep_end}        = $poly_len;
		$$pep{pep_sequence}   = substr( $poly_protein, $poly_pos );
		$$pep{pep_length}     = length( $$pep{pep_sequence} );
		if ( $$pep{pep_sequence} =~ /\*^/ ) { $$pep{pep_length}-- }
		push @matpeps, $pep;
	}

	return \@matpeps;
}

sub compare_cleavage_sites {
	my ( $a, $b ) = @_;

	my ( undef, $offseta, undef  ) = split / *\| */, $a;
	my ( undef, $offsetb, undef  ) = split / *\| */, $b;

	return $offseta <=> $offsetb;
}

sub merge_partial_genes {
	my ( $genome, $genes ) = @_;

	my @merged;
	my $verbose = get_parameter( "verbose" );
	#$verbose = 1;

	# find genes fragmented into partial genes
	my %gene_syms;
	for my $partial ( @$genes ) {
		my $gene_sym = "$$partial{gene_name}|$$partial{orientation}";
		if ( exists $gene_syms{$gene_sym} ) {
			if ( $$partial{is_mutation} && ! $gene_syms{$gene_sym}{is_mutation} ) { next }
			$gene_syms{$gene_sym}{is_mutation} = $$partial{is_mutation};
			push @{$gene_syms{$gene_sym}{partials}}, $partial;
		}
		else {
			my @partials = ( $partial );
			$gene_syms{$gene_sym}{is_mutation} = $$partial{is_mutation};
			$gene_syms{$gene_sym}{partials} = \@partials;
		}
	}

	# merge partials into single gene
	for my $gene_sym ( keys %gene_syms ) {

		# mutation or unfragmented gene, use as is
		my @partials = sort { compare_gene_positions( $a, $b ) } @{$gene_syms{$gene_sym}{partials}};

		# nothing to merge, add partials "as is"
		if ( @partials == 1 || $gene_syms{$gene_sym}{is_mutation} ) {
			push @merged, @partials;
		}

		# gene fragmented into partials, merge all pieces into single gene
		else {

			# initialize gene from first partial instance
			if ( $verbose ) {
				print "\nMERGE partials for $gene_sym\n";
				print_genehits( @partials );
			}
			my $gene = shift @partials;
			my %refs;
			$refs{$$gene{ref_id}} = 1;

			# merge remaining partials into gene
			while ( @partials ) {
				my $partial = shift @partials;
				$refs{$$partial{ref_id}} = 1;

				# merge exons, cdna, and protein
				my $exons = $$gene{exons};
				my $next_exons = $$partial{exons};
				my $cdna_adjustment = length( $$gene{cdna} );

				# fuse last/first?
				if ( ! allow_splicing($gene) ) {
					my $gapstart = $$gene{stop_site} + $$gene{orientation};
					my $gapend = $$partial{start_codon} - $$gene{orientation};
					my $tweendna = subsequence( $$genome{sequence}, $gapstart, $gapend );
					my $tweenlen = length( $tweendna );
					my $tweenpep = DNA2AA( $tweendna );

					# check for fusion errors
					my $skip_fusion = 0;
					if ( $tweenlen % 3 != 0 ) {
						if ( allow_ribosomal_slippage($gene) ) {
							$skip_fusion = 1;
						}
						else {
							my $msg = "gap induced frameshift between "
								. minval( $gapstart, $gapend ) . " and " . maxval( $gapstart, $gapend );
							if ( ! $$gene{is_mutation} ) {
								$$gene{is_mutation} = 1;
								$$gene{warning} = $msg;
							}
							else {
								if ( $$gene{warning} !~ /$msg/ ) {
									$$gene{warning} .= ", $msg";
								}
							}
						}
					}

					if ( ! $skip_fusion ) {
						my $stop = index( $tweenpep, "*" );
						if ( $stop >= 0 ) {
							my $stopstart = $gapstart + 3 * $stop * $$gene{orientation};
							my $stopend = $stopstart + 2 * $$gene{orientation};
							my $msg = "ORF interrupted by stop codon at $stopstart..$stopend";
							if ( ! $$gene{is_mutation} ) {
								$$gene{is_mutation} = 1;
								$$gene{warning} = $msg;
							}
							else {
								if ( $$gene{warning} !~ /$msg/ ) {
									$$gene{warning} .= ", $msg";
								}
							}
						}
						if ( $$gene{protein} =~ /\*$/ ) {
							my $stopend = $$gene{stop_site};
							my $stopstart = $stopend - 2 * $$gene{orientation};
							my $msg = "ORF interrupted by stop codon at $stopstart..$stopend";
							if ( ! $$gene{is_mutation} ) {
								$$gene{is_mutation} = 1;
								$$gene{warning} = $msg;
							}
							else {
								if ( $$gene{warning} !~ /$msg/ ) {
									$$gene{warning} .= ", $msg";
								}
							}
						}

						# fuse last exon of current gene with first exon of this partial
						$$gene{cdna} .= $tweendna;
						$cdna_adjustment += $tweenlen;
						my $next = shift @$next_exons;
						my $last = $$exons[@$exons-1];
						$$last{dna_end} = $$next{dna_end};
						$$last{cdna_end} = $$next{cdna_end} + $cdna_adjustment;
						$$gene{protein} .= $tweenpep;
					}
				}

				# merge cdna, exons, and proteins
				$$gene{cdna} .= $$partial{cdna};
				for my $exon ( @$next_exons ) {
					$$exon{cdna_start} += $cdna_adjustment;
					$$exon{cdna_end} += $cdna_adjustment;
				}
				push @$exons, @$next_exons;
				$$gene{protein} .= $$partial{protein};
				$$gene{protein_length} = length( $$gene{protein} );
				if ( $$gene{protein} =~ /\*$/ ) { $$gene{protein_length}-- }

				# extend end of gene
				$$gene{stop_site} = $$partial{stop_site};
				$$gene{stop_truncated} = $$partial{stop_truncated};

#				if ( $verbose && @partials ) {
#					print "INTERMEDIATE\n";
#					print_genehits( $gene );
#					print "gene=$$gene{protein}\n";
#				}
			}

			# update reference sequence
			$$gene{match_quality} = 0;
			for my $ref_id ( keys %refs ) {
				my %tmpgene = %$gene;
				$tmpgene{ref_id} = $ref_id;
				$tmpgene{gene_definition} = get_reference_definition( $ref_id );
				score_merged_partials( $genome, \%tmpgene );
				if ( $tmpgene{match_quality} > $$gene{match_quality} ) { %$gene = %tmpgene }
			}
			finalize_merged_partials( $genome, $gene );
			validate_reference_alignment( $genome, $gene );

			# save merged gene
			push @merged, $gene;
			if ( $verbose ) {
				print "\nMERGED\n";;
				print_genehits( $gene );
#				print "gene=$$gene{protein}\n";
			}
		}
	}

	@$genes = sort { compare_gene_ids( $a, $b ) } @merged;
}

####################################
# pairwise alignments
# align peptide to reference
sub align_pep_to_reference {
	my ( $pepid, $pepseq, $refid, $refseq ) = @_;

	my $vigorspace = get_parameter("vigorspace");

	my $tmppepid = $pepid;
	$tmppepid =~ s/:/_/g;
	my $tmprefid = $refid;
	$tmprefid =~ s/:/_/g;

	open( TMP, ">$vigorspace/tmp_align_fasta" )
	  || die "Could not write to the tmp_align_fasta file\.\n";
	print TMP ">$tmppepid\n$pepseq\n>$tmprefid\n$refseq\n";
	close TMP;

	my $alignment;
	my $cmd = "clustalw -align -infile=$vigorspace/tmp_align_fasta &> $vigorspace/tmp_align_fasta.log";
	if ( !system $cmd ) {
        sleep(2);
		$alignment = `cat $vigorspace/tmp_align_fasta.aln`;
		if ( $tmppepid ne $pepid ) {
			$tmppepid =~ s/([\.\/\|\\])/\\$1/g;
			$alignment =~ s/$tmppepid/$pepid/g;
		}
		if ( $tmppepid ne $pepid ) {
			$tmprefid =~ s/(\W)/\\$1/g;
			$alignment =~ s/$tmprefid/$refid/g;
		}
		$alignment =~ s/[\n\r\s]+$//;
#print "alignment $pepid vs $refid=\n$alignment\n";
	}
	else {
		$alignment = undef;
	}

	unlink "$vigorspace/tmp_align_fasta";
	unlink "$vigorspace/tmp_align_fasta.aln";
	unlink "$vigorspace/tmp_align_fasta.log";

	return $alignment;
}

# align mutation to reference
sub bl2seq_alignment {
	my ( $genome, $gene, $bl2seqopts ) = @_;
	my $vigorspace = get_parameter("vigorspace");

	# get sequence for extended gene span
	my $genome_sequence = $$genome{sequence};
	my $extend_span = 50;
	my $detailed_exceptions = get_parameter("detailed_exceptions");

	my $gene_span = "";

	my $first = minval( $$gene{start_codon}, $$gene{stop_site} );
	my $left = $first - $extend_span;
	if ( $left < 1 ) {
		$left = 1;
	}
	my $right = $left + $extend_span + abs( $$gene{stop_site} - $$gene{start_codon} ) + $extend_span ;
	if ( $right > length($genome_sequence) ) {
		$left -= $right - length($genome_sequence);
		if ( $left < 1 ) { $left = 1 }
		$right = length($genome_sequence);
	}

	my $end = $left - 1;
	for my $exon( sort { $$a{dna_start} <=> $$b{dna_start} } @{$$gene{exons}} ) {
		my $start = minval( $$exon{dna_start}, $$exon{dna_end} );
		if ( $start > $end+1 ) {
			$gene_span .= lc subsequence( $genome_sequence, $end+1, $start - 1 );
		}
		elsif ( $start <= $end ) {
			$start = $end+1;
		}
		$end = maxval( $$exon{dna_start}, $$exon{dna_end} );
		if ( $end >= $start ) { $gene_span .= uc subsequence( $genome_sequence, $start, $end ) }
	}

	if ( $right > $end ) {
		$gene_span .= lc subsequence( $genome_sequence, $end + 1, $right );
	}

	if ( $$gene{orientation} == -1 ) { $gene_span = reverse_complement( $gene_span ) }

	$gene_span =~ s/(.{60})/$1\n/g;
	if ( substr( $gene_span, length($gene_span) - 1, 1 ) ne "\n" ) {
		$gene_span .= "\n";
	}

	my $dna_alignment = ">$$gene{gene_id}\t$$gene{start_codon}\t$$gene{stop_site}\n";
	$dna_alignment .= $gene_span . "\n";

	# use bl2seq to align extended gene span to reference sequence
	unlink "$vigorspace/mutantspan.fasta";
	open( TMP, ">$vigorspace/mutantspan.fasta" )
	  || die "could not write mutantspan.fasta to disk";
	print TMP ">$$gene{gene_id}\n$gene_span";
	close TMP;

	unlink "$vigorspace/mutantref.fasta";
	my $ref = get_reference_seq( $$gene{ref_id} );
	my $ref_sequence = $$ref{sequence};
	open( TMP, ">$vigorspace/mutantref.fasta" )
	  || die "could not write mutantref.fasta to disk";
	print TMP ">$$ref{defline}\n$ref_sequence";
	close TMP;

	my $cmd;
	$cmd = "bl2seq $bl2seqopts -i $vigorspace/mutantspan.fasta -j $vigorspace/mutantref.fasta";
	my @tmp = split /\n/, `$cmd`;
	#print " ref=" . `cat $vigorspace/mutantref.fasta` . "\n";
	#print "span=" . `cat $vigorspace/mutantspan.fasta` . "\n";
	#print "\nBL2SEQ ------------------------------\ncmd=$cmd\n" . join( "\n", @tmp ) . "\n-----------------------------------\n";
	# we need to adjust the coordinates back to the full genomic sequence

	# calculate the largest possible coordinate so we can size the margin appropriately
	my $coordsz = $$gene{start_codon};
	if ( $$gene{stop_site} > $coordsz ) { $coordsz = $$gene{stop_site} }
	if ( length( $reference_seqs{ $$gene{ref_id} }{sequence} ) > $coordsz )
	{
		$coordsz = length( $reference_seqs{ $$gene{ref_id} }{sequence} );
	}
	$coordsz = length($coordsz);

	# load the alignment results
	# discard junk at start of alignment output
	while ( @tmp && $tmp[0] !~ /^  *Length = / ) {
#if ( $$gene{is_mutation} ) { print "skip $tmp[0]\n" }
		shift @tmp;
	}
	#shift @tmp;

	# calculate offset adjustment for the coordinates
	my $offset = $left - 1;
	if ( $$gene{orientation} == -1 ) { $offset = $right+1 }

	# reformat the alignment with adjusted coordinates
	my $oldmargin;
	my $newmargin;
	for my $line (@tmp) {
#if ( $$gene{is_mutation} ) { print "line $line\n" }

		# discard everything from "Lambda" on
		if ( $line =~ /Lambda  *K  *H/ ) { last }

		# adjust query alignment string
		if ( $line =~ /(Query: +([0-9]+) +)([^ ]+) +([0-9]+) *$/ ) {
			my ( $prefix, $start, $seq, $end ) = ( $1, $2, $3, $4 );
			$oldmargin = length($prefix);
			if ( $$gene{orientation} == 1 ) {
				$start = $offset + $start;
				$end   = $offset + $end;
			}
			else {
				$start = $offset - $start;
				$end   = $offset - $end;
			}
			$prefix =
			  "Query: " . substr( "$start            ", 0, $coordsz ) . " ";
			if ( !defined $newmargin ) { $newmargin = length($prefix) }
			$line = $prefix . $seq . " " . $end;
		}

		# adjust subject alignment string
		elsif ( $line =~ /(Sbjct: +([0-9]+) +)([^ ]+) +([0-9]+) *$/ ) {
			my ( $prefix, $start, $seq, $end ) = ( $1, $2, $3, $4 );
			$prefix =
			  "Sbjct: " . substr( "$start            ", 0, $coordsz ) . " ";
			$line = $prefix . $seq . " " . $end;
		}

		# adjust middle alignment string (comparision of query to subject)
		elsif ( defined $newmargin && $line !~ /^\s*$/ ) {
			if ( $newmargin < $oldmargin ) {
				$line = substr( $line, $oldmargin - $newmargin );
			}
			elsif ( $newmargin > $oldmargin ) {
				$line =
				  substr( "            ", 0, $newmargin - $oldmargin ) . $line;
			}
		}

		# add line to final result
		if ( $line =~ /^ ((Sbjct|Query)\: [0-9]+)([^ 0-9].*$)/ ) { $line = "$1 $3"}
		$dna_alignment .= "$line\n";
	}

	return $dna_alignment;
}

########################################
# reports

# write standard TBL file
sub std_tbl_report {
	my ( $TBL, $genome, $genes ) = @_;

	if ( ! defined $TBL ) { return }

	if ( ! defined $genes ) { return }

	my $feature_count = 0;
	for my $gene ( @$genes ) {
		if ( $$gene{is_mutation} ) {
			my $location = get_gene_location( $gene );
			my @tmp = split /\.\./, $location;
			my $start = shift @tmp;
			my $end = pop @tmp;
			my $subject_coverage = ( $end - $start + 1 ) / length( $$genome{sequence} );
			if ( $$gene{pct_refcoverage} >= 95.0 || $$gene{pct_refsimilarity} >= 95.0 ) {
				$feature_count++;
			}
		}
		else {
			$feature_count++;
		}
	}
	if ( $feature_count == 0 ) { return }

	print $TBL ">Features $$genome{id}\n";

	for my $gene ( @$genes ) {
		my $location = get_gene_location( $gene );
		if ( $$gene{is_mutation} ) {
			my @tmp = split /\.\./, $location;
			my $start = shift @tmp;
			my $end = pop @tmp;
			my $subject_coverage = ( $end - $start + 1 ) / length( $$genome{sequence} );
			if ( $$gene{pct_refcoverage} >= 95.0 || $$gene{pct_refsimilarity} >= 95.0 ) {
				print $TBL "$start\t$end\tfeature\n";
				print $TBL "\t\t\tnote\tsimilar to $$gene{gene_definition}, $$gene{warning}\n";
			}
		}
		else {
			my @tmp = split /\.\./, $location;
			my $start = shift @tmp;
			my $end = pop @tmp;
			print $TBL "$start\t$end\tgene\n";
			print $TBL "\t\t\tgene\t$$gene{gene_name}\n";

			my @exons = split /,/, $location;
			for my $i ( 0 .. @exons-1 ) {
				my $exon = $exons[$i];
				my ( $exstart, $exend ) = split /\.\./, $exon;
				if ( $i == 0 ) {
					print $TBL "$exstart\t$exend\tCDS\n";
				}
				else {
					print $TBL "$exstart\t$exend\n";
				}
			}
			if ( exists $$gene{translation_exception} ) {
				print $TBL "\t\t\ttransl_except\t(pos:$$gene{translation_exception}{begin}..$$gene{translation_exception}{end},aa:$$gene{translation_exception}{aa})\n";
			}
			if ( allow_ribosomal_slippage($gene) && @exons==2 ) {
				print $TBL "\t\t\tribosomal_slippage\n";
			}
			print $TBL "\t\t\tprotein_id\t$$gene{gene_id}\n";
			print $TBL "\t\t\tgene\t$$gene{gene_name}\n";
			print $TBL "\t\t\tproduct\t$$gene{gene_definition}\n";
			if ( $$gene{warning} ) {
				print $TBL "\t\t\tnote\t$$gene{warning}\n";
			}
		}

	}
	for my $gene ( @$genes ) {
		if ( defined $$gene{mature_peps} ) {
			print $TBL "\n>Features $$gene{gene_id}\n";
			for my $pep ( @{ $$gene{mature_peps} } ) {
				my $start = $$pep{pep_start};
				if ( $$gene{start_truncated} && $start <= 1 ) {
					$start = "<$start";
				}
				my $end = $$pep{pep_end};
				if ( $$gene{stop_truncated} && $end >= $$gene{protein_length} ) {
					$end = ">$end";
				}
				print $TBL "$start\t$end\tmat_peptide\n";
				print $TBL "\t\t\tgene\t$$gene{gene_name}\n";
				print $TBL "\t\t\tproduct\t$$pep{pep_definition}\n";
			}
		}
	}
}

# write cds fasta
sub std_cds_report {
	my ( $CDS, $genome, $genes ) = @_;

	if ( ! defined $CDS ) { return }

	for my $gene ( @$genes ) {
		if ( ! $$gene{is_mutation} ) {

			my $location = get_gene_location( $gene );
			my $cdna = $$gene{cdna};

			my $defline = ">$$gene{gene_id} location=$location";
			if ( exists $$gene{translation_exception} ) {
				$defline .= " translation_exception=\"$$gene{translation_exception}{begin}..$$gene{translation_exception}{end}:$$gene{translation_exception}{aa}\"";
			}
			$defline .= " gene=\"$$gene{gene_name}\" product=\"$$gene{gene_definition}\"\n";
			print $CDS $defline;

			$cdna =~ s/(.{60})/$1\n/g;
			if ( substr( $cdna, length($cdna) - 1, 1 ) ne "\n" ) { $cdna .= "\n" }

			print $CDS  $cdna . "\n";

			if ( defined $$gene{mature_peps} ) {
				for my $pep ( @{ $$gene{mature_peps} } ) {

					my $location = get_matpep_location( $gene, $pep );

					my @tmp_description;
					if ( defined $$pep{pep_name} ) {
						push @tmp_description, $$pep{pep_name};
					}
					if ( defined $$pep{pep_definition} ) {
						push @tmp_description, $$pep{pep_definition};
					}
					my $pep_description = join( " ", @tmp_description );
					if ( $pep_description !~ /mat[ure]* *pep[tide]*/i ) {
						$pep_description = "mature peptide, $pep_description";
					}

					my $defline = ">$$pep{pep_id} location=$location gene=\"$$gene{gene_name}\" product=\"$pep_description\"\n";

					print $CDS $defline;

					my $cdna = "";
					for my $exon ( split /,/, $location ) {
						$exon =~ s/[><]//g;
						my( $start, $end ) = split /\.\./, $exon;
						$cdna = subsequence( $$genome{sequence}, $start, $end );
					}
					$cdna =~ s/(.{60})/$1\n/g;
					if ( substr( $cdna, length($cdna) - 1, 1 ) ne "\n" ) {
						$cdna .= "\n";
					}

					print $CDS $cdna . "\n";
				}
			}
		}
	}
}

# write peptide fasta
sub std_pep_report {
	my ( $PEP, $genome, $genes ) = @_;

	if ( ! defined $PEP ) { return }

	for my $gene ( @$genes ) {
		if ( ! $$gene{is_mutation} ) {

			my $location = get_gene_location( $gene );
			my $protein = $$gene{protein};
			my $numaa = length( $protein );
			if ( substr( $protein, $numaa-1, 1 ) eq "*" ) {
				$numaa--;
			}

			my $defline = ">$$gene{gene_id} location=$location";
			if ( exists $$gene{translation_exception} ) {
				$defline .= " translation_exception=\"$$gene{translation_exception}{begin}..$$gene{translation_exception}{end}:$$gene{translation_exception}{aa}\"";
			}
			$defline .= " length=$numaa gene=\"$$gene{gene_name}\" product=\"$$gene{gene_definition}\"\n";

			print $PEP $defline;

			$protein =~ s/(.{60})/$1\n/g;
			if ( substr( $protein, length($protein) - 1, 1 ) ne "\n" ) { $protein .= "\n" }

			print $PEP $protein . "\n";

			if ( defined $$gene{mature_peps} ) {
				for my $pep ( @{ $$gene{mature_peps} } ) {

					my $protein = $$pep{pep_sequence};
					my $numaa = length( $protein );
					if ( substr( $protein, $numaa-1, 1 ) eq "*" ) {
						$numaa--;
					}

					my $location = get_matpep_location( $gene, $pep );

					my @tmp_description;
					if ( defined $$pep{pep_name} ) {
						push @tmp_description, $$pep{pep_name};
					}
					if ( defined $$pep{pep_definition} ) {
						push @tmp_description, $$pep{pep_definition};
					}
					my $pep_description = join( " ", @tmp_description );

					print $PEP ">$$pep{pep_id} location=$location length=$numaa gene=\"$$gene{gene_name}\" product=\"$pep_description\"\n";

					$protein =~ s/(.{60})/$1\n/g;
					if ( substr( $protein, length($protein) - 1, 1 ) ne "\n" ) {
						$protein .= "\n";
					}

					print $PEP $protein . "\n";
				}
			}
		}
	}
}

# write standard ALIGN report
sub std_align_report {
	my ( $ALIGN, $FS, $genome, $genes ) = @_;

	if ( ! defined $ALIGN && ! defined $FS ) { return }

	for my $gene ( @$genes ) {
		if ( defined $ALIGN ) { gene_align_report( $ALIGN, $gene ) }
		if ( defined $FS &&  defined $$gene{warning} && $$gene{warning} =~ /frameshift/i ) {
			gene_align_report( $FS, $gene );
		}
	}
	return;
}

sub gene_align_report {
	my ( $OUT, $gene ) = @_;

	if ( ! defined $OUT ) { return }

	for my $i ( 1..24 ) { print $OUT "=====" }
	print $OUT "\n\n";

	print $OUT format_genehits( $gene ) . "\n";
	if ( defined $$gene{pep_alignment} ) {
		print $OUT "$$gene{pep_alignment}\n";
		print $OUT "\n";
	}

	if ( defined $$gene{ref_alignment} ) {
		print $OUT "$$gene{ref_alignment}\n";
		print $OUT "\n";
	}

	if ( defined $$gene{mp_alignment} ) {
		print $OUT "$$gene{mp_alignment}\n";
		print $OUT "\n";
	}
}

sub print_blasthits {
	my ( @blasthits ) = @_;
	my $blanks40 = "                                        ";
	my $header =
	    substr( "subject_id$blanks40", 0, 17 ) . " "
	  . substr( "ori$blanks40", 0, 3 ) . " "

	  #	  . substr( "subject_id$blanks40", 0, 14 ) . " "
	  . substr( "subject_range$blanks40", 0, 14 ) . " "
	  . substr( "query_range$blanks40", 0, 14 ) . " "
	  . substr( "qf$blanks40",   0, 2 ) . " "

	  #	  . substr( "ppenalty$blanks40", 0, 8 ) . " "
	  #	  . substr( "frame$blanks40", 0, 5 ) . " "
	  . substr( "hsp",             0, 3 ) . " "
	  . substr( "alglen$blanks40", 0, 6 ) . " "
	  . substr( "lorf$blanks40",   0, 4 ) . " "
	  . substr( "stp$blanks40",    0, 3 ) . " "
	  . substr( "fs$blanks40",     0, 2 ) . " "
	  . substr( "ident$blanks40",  0, 5 ) . " "
	  . substr( "simil$blanks40",  0, 5 ) . " "
	  . substr( "%id$blanks40",    0, 6 ) . " "
	  . substr( "%slen$blanks40",  0, 6 ) . " "
	  . substr( "%vsim$blanks40",  0, 6 ) . " "
	  . substr( "vwgt$blanks40",   0, 6 ) . " "
	  . substr( "evalue$blanks40", 0, 8 ) . " "
	  . "subject definition";
	print "$header\n";

	my $oldsid = "";
	for my $blasthit (@blasthits) {

		my @hits = ( $blasthit );
		if ( exists $$blasthit{permutations} ) {
			push @hits, @{$$blasthit{permutations}};
		}

		for my $hit ( @hits ) {
			my $subject_id = substr( "$$hit{subject_id}$blanks40", 0, 17 );

			my $query_range =
			  substr( "$$hit{query_left}-$$hit{query_right}$blanks40", 0, 14 );

			my $query_frame =
			  substr( "$$hit{query_frame}$blanks40", 0, 2 );

			my $subject_range =
			  substr( "$$hit{subject_left}-$$hit{subject_right}$blanks40", 0, 14 );

			my $ori = substr( "$$hit{orientation}$blanks40", 0, 3 );

			#		my $frame = substr( "$$hit{query_frame}$blanks40", 0, 5 );

			my $hsps;
			if ( defined $$hit{hsps} ) {
				$hsps = @{ $$hit{hsps} };
			}
			else {
				$hsps = 0;
			}
			$hsps = substr( "$hsps$blanks40", 0, 3 );

			my $alignment_length =
			  substr( "$$hit{alignment_length}$blanks40", 0, 6 );

			my $lorf = substr( "$$hit{longest_orf}$blanks40", 0, 4 );

			my $stops = substr( "$$hit{query_stops}$blanks40", 0, 3 );

			my $frameshifts = substr( "$$hit{has_frameshift}$blanks40", 0, 2 );

			my $num_identical = substr( "$$hit{num_identical}$blanks40", 0, 5 );

			my $num_similar = substr( "$$hit{num_similar}$blanks40", 0, 5 );

			my $pct_identity = $$hit{pct_identity};
			if ( index( $pct_identity, "." ) < 0 ) { $pct_identity .= ".0" }
			$pct_identity = substr( "$pct_identity$blanks40", 0, 6 );

			my $pct_scoverage = $$hit{pct_scoverage};
			if ( index( $pct_scoverage, "." ) < 0 ) { $pct_scoverage .= ".0" }
			$pct_scoverage = substr( "$pct_scoverage$blanks40", 0, 6 );

			my $pct_vsimilarity = int( 10.0 * $$hit{vigor_pctsimilarity} ) / 10.0;
			if ( index( $pct_vsimilarity, "." ) < 0 ) { $pct_vsimilarity .= ".0" }
			$pct_vsimilarity = substr( "$pct_vsimilarity$blanks40", 0, 6 );

			my $evalue = substr( format_evalue( $$hit{evalue} ) . $blanks40, 0, 8 );

			my $vweight = int( 10.0 * $$hit{vigor_matchwgt} ) / 10.0;
			if ( index( $vweight, "." ) < 0 ) { $vweight .= ".0" }
			$vweight = substr( "$vweight$blanks40", 0, 6 );

			my $definition = get_reference_definition( $$hit{subject_id} );
			if ( ! defined $definition ) {
				$definition = $$hit{subject_definition};
				$definition =~ s/^[^ ]* *//;
			}

			my $line = join(
				" ",
				( $subject_id, $ori, $subject_range, $query_range, $query_frame, $hsps, $alignment_length,
				  $lorf, $stops, $frameshifts, $num_identical, $num_similar, $pct_identity,
				  $pct_scoverage, $pct_vsimilarity, $vweight, $evalue, $definition
				)
			);
			print "$line\n";
		}
	}
	return;
}

sub stats_report{
	my ( $STATS, $genome, $genelist ) = @_;

	my $frameshifts = 0;
	my $mutations = 0;
	my $genes = 0;
	my $cdsbases = 0;
	my $refbases = 0;
	my $pepbases = 0;
	my $coveredaa = 0;
	my $identicalaa = 0;
	my $similaraa = 0;

	for my $gene ( @$genelist ) {
		if ( defined $$gene{warning} && $$gene{warning} =~ /frameshift/i ) { $frameshifts++; }
		if ( $$gene{is_mutation} ) {
			$mutations++;
		}
		else {
			$genes++;
			$cdsbases += length( $$gene{cdna} );
			$refbases += $$gene{ref_length};
			$pepbases += $$gene{protein_length};
			$coveredaa += int( $$gene{pct_refcoverage}/100.0 * $$gene{ref_length} + 0.5 );
			$identicalaa += int( $$gene{pct_refidentity}/100.0 * $$gene{ref_length} + 0.5 );
			$similaraa += int( $$gene{pct_refsimilarity}/100.0 * $$gene{ref_length} + 0.5 );
		}
	}
	my $coverage = "-";
	my $similarity = "-";
	my $identity = "-";
	if ( $refbases ) {
		$coverage = int( 1000.0 * $coveredaa / $refbases + 0.5 ) / 10.0;
		if ( index( $coverage, ".") < 0 ) { $coverage .= ".0" }
		$identity = int( 1000.0 * $identicalaa / $refbases + 0.5 ) / 10.0;
		if ( index( $identity, ".") < 0 ) { $identity .= ".0" }
		$similarity = int( 1000.0 * $similaraa / $refbases + 0.5 ) / 10.0;
		if ( index( $similarity, ".") < 0 ) { $similarity .= ".0" }
	}

	print join( "  ",
			( rpad( $$genome{id}, 20),
				lpad( length( $$genome{sequence} ), 7 ),
				lpad( $genes, length("Genes") ),
				lpad( $mutations, length("Errors") ),
				lpad( $frameshifts, length("Frameshifts") ),
				lpad( $cdsbases, length("CDS Bases") ),
				lpad( $pepbases, length("Peptide Bases") ),
				lpad( $coverage, length("%Ref Coverage") ),
				lpad( $identity, length("%Ref Identity") ),
				lpad( $similarity, length("%Ref Similarity") ) ) )
		. "\n";

	if ( defined $STATS ) {
		print $STATS join( "\t",
				( $$genome{id}, length( $$genome{sequence} ), $genes, $mutations, $frameshifts, $cdsbases, $pepbases, $coverage, $identity, $similarity ) )
			. "\n";
	}
}

sub stats_header{
	my ( $STATS, $STDOUT ) = @_;

	if ( defined $STDOUT ) {
		print $STDOUT join( "  ",
				( rpad( "Sequence", 20),
					lpad( "Length", 7),
					lpad( "Genes", length("Genes") ),
					lpad( "Errors", length("Errors") ),
					lpad( "Frameshifts", length("Frameshifts") ),
					lpad( "CDS Bases", length("CDS Bases") ),
					lpad( "Peptide Bases", length("Peptide Bases") ),
					lpad( "%Ref Coverage", length("%Ref Coverage") ),
					lpad( "%Ref Identity", length("%Ref Identity") ),
					lpad( "%Ref Similarity", length("%Ref Similarity") ) ) )
			. "\n";
	}

	if ( defined $STATS ) {
		print $STATS join( "\t",
				( "Sequence", "Length", "Genes", "Errors", "Frameshifts", "CDS Bases", "Peptide Bases", "%Ref Coverage", "%Ref Identity", "%Ref Similarity" ) )
			. "\n";
	}
}

sub print_genehits {
	my ( @genes ) = @_;
	print format_genehits( @genes );
}

sub format_genehits {
	my ( @genes ) = @_;

	my $blanks50 = "                                                  ";
	my $header =
	  substr( "gene_id$blanks50", 0, 18 ) . " "
	  . "Err "
	  . substr( "%cov$blanks50",       0, 6 ) . " "
	  . substr( "%id$blanks50",        0, 6 ) . " "
	  . substr( "%sim$blanks50",       0, 6 ) . " "
	  . "Ex "
	  . substr( "start..stop$blanks50", 0, 20 ) . " "
	  . substr( "%t5$blanks50",        0, 5 ) . " "
	  . substr( "%t3$blanks50",        0, 5 ) . " "
	  . substr( "ref_id$blanks50",     0, 12 ) . " "
	  . "definition";

	my $hits = "$header\n";

	for my $gene (@genes) {
		my $gene_id = substr( "$$gene{gene_id}$blanks50", 0, 18 );

		my $mutation = " $$gene{is_mutation} ";

		my $cov = "      ";
		if ( defined $$gene{pct_refcoverage} ) {
			$cov = substr( "$$gene{pct_refcoverage}      ", 0, 6 );
		}

		my $ident = "      ";
		if ( defined $$gene{pct_refidentity} ) {
			$ident = substr( "$$gene{pct_refidentity}      ", 0, 6 );
		}

		my $sim = "      ";
		if ( defined $$gene{pct_refsimilarity} ) {
			$sim = substr( "$$gene{pct_refsimilarity}      ", 0, 6 );
		}

		my $num_exons = @{ $$gene{exons} };
		$num_exons = substr( "$num_exons   ", 0, 2 );

		my $query_range = get_gene_location( $gene );
		my @ranges = split /,/, $query_range;
		$query_range = shift @ranges;
		$query_range = substr( "$query_range$blanks50", 0, 20 );

		my $t5 = "     ";
		if ( defined $$gene{pct_reftrunc5} ) {
			$t5 = substr( "$$gene{pct_reftrunc5}      ", 0, 5 );
		}

		my $t3 = "     ";
		if ( defined $$gene{pct_reftrunc3} ) {
			$t3 = substr( "$$gene{pct_reftrunc3}      ", 0, 5 );
		}

		my $subject_id = $$gene{ref_id};
		$subject_id = substr( "$subject_id$blanks50", 0, 12 );

		my $definition = $$gene{gene_name} . " | " . $$gene{gene_definition};
		my $line = join( " ",
			( $gene_id, $mutation, $cov, $ident, $sim,
			$num_exons, $query_range, $t5, $t3,
			$subject_id, $definition ) );
		$hits .= "$line\n";

		for my $range ( @ranges ) {
			$hits .=
				substr( $blanks50, 0, 18 ) . " "			# gene_id
				. "    "									# Err
				. substr( $blanks50, 0, 6 ) . " "			# %cov
				. substr( $blanks50, 0, 6 ) . " "			# %id
				. substr( $blanks50, 0, 6 ) . " "			# %sim
				. "   "										# Ex
				. $range									# start..stop
				. "\n";
		}
		if ( $$gene{warning} ) {
			$hits .= substr( $blanks50, 0, 20 ) . $$gene{warning} . "\n";
		}

		if ( defined $$gene{mature_peps} ) {
			for my $pep ( @{$$gene{mature_peps}} ) {
				my $location = get_matpep_location( $gene, $pep );
				my @ranges = split /,/, $location;

				my @tmp_description;
				if ( defined $$pep{pep_name} ) {
					push @tmp_description, $$pep{pep_name};
				}
				if ( defined $$pep{pep_definition} ) {
					push @tmp_description, $$pep{pep_definition};
				}
				my $pep_description = join( " ", @tmp_description );

				my $range = shift @ranges;

				if ( defined $$pep{ref_id} ) {
					my $pep_id = substr( "$$pep{pep_id}$blanks50", 0, 20 );

					$cov = "      ";
					if ( defined $$pep{pct_refcoverage} ) {
						$cov = substr( "$$pep{pct_refcoverage}      ", 0, 6 );
					}

					$ident = "      ";
					if ( defined $$pep{pct_refidentity} ) {
						$ident = substr( "$$pep{pct_refidentity}      ", 0, 6 );
					}

					$sim = "      ";
					if ( defined $$pep{pct_refsimilarity} ) {
						$sim = substr( "$$pep{pct_refsimilarity}      ", 0, 6 );
					}

					$num_exons = "  ";

					$query_range = substr( "$range$blanks50", 0, 20 );

					$t5 = "     ";
					if ( defined $$pep{pct_reftrunc5} ) {
						$t5 = substr( "$$pep{pct_reftrunc5}      ", 0, 5 );
					}

					$t3 = "     ";
					if ( defined $$pep{pct_reftrunc3} ) {
						$t3 = substr( "$$pep{pct_reftrunc3}      ", 0, 5 );
					}

					$subject_id = substr( "$$pep{ref_id}$blanks50", 0, 12 );;

					my $line = join( " ", ( $pep_id, " ",  $cov, $ident, $sim,
						$num_exons, $query_range, $t5, $t3, $subject_id, $pep_description ) );

					$hits .= "$line\n";
				}
				else {
					$hits .=
						substr( "$$pep{pep_id}$blanks50", 0, 47 )	# pep_id
						. substr( "$range$blanks50", 0, 46 )		# start..stop
	  					. $pep_description							# definition
						. "\n";
				}

				for $range ( @ranges ) {
					$hits .=
						substr( $blanks50, 0, 47 )				# pep_id
						. $range								# start..stop
						. "\n";
				}

			}
		}
	}

	return $hits;
}

sub print_no_genes {
	my ( $genomic_seq, $ALIGN ) = @_;

	for my $i ( 1..24 ) { print $ALIGN "=====" }

	print $ALIGN ">$$genomic_seq{id}\n\n"
	  . "ERROR MESSAGES\n\n"
	  . "Wrong Virus Type Was Selected Or Sequence Is Too Short Or Sequences Were Too Divergent.\n\n"
	  . "########\n";

	print ">$$genomic_seq{id}\n\n"
	  . "ERROR MESSAGES\n\n"
	  . "Wrong Virus Type Was Selected Or Sequence Is Too Short Or Sequences Were Too Divergent.\n\n"
	  . "########\n";
}

sub get_blasthit_location {
	my ( $hit ) = @_;

	my @hsps;
	if ( ! exists $$hit{hsps} ) {
		@hsps = ( $hit );
	}
	else {
		@hsps = @{$$hit{hsps}}
	}

	my $location = "";
	for my $hsp ( sort { $$a{orientation} * $$a{query_left} <=> $$b{orientation} * $$b{query_left} } @hsps ) {
		if ( length( $location ) ) { $location .= "," }
		if ( $$hsp{orientation} == 1 ) {
			$location .= "$$hsp{query_left}..$$hsp{query_right}";
		}
		else {
			$location .= "$$hsp{query_right}..$$hsp{query_left}";
		}
	}

	return $location;
}

sub get_gene_location {
	my ( $gene ) = @_;

	my $dbg = 0;

	if ( $$gene{is_mutation} ) {
		return "$$gene{start_codon}..$$gene{stop_site}";
	}

	my $location =
		cdna_coords_to_dna( $$gene{orientation}, $$gene{exons}, 1, length( $$gene{cdna} ), 1 );

	if ( $dbg ) {
		print "start $$gene{start_codon}  stop $$gene{stop_site}";
		for my $exon ( @{$$gene{exons}} ) {
			print "  exon $$exon{dna_start}-$$exon{dna_end}"
		}
		print "  length " . length( $$gene{cdna} ) . "  location $location\n";
	}
	my $allow_partial = get_parameter( "allow_partial_genes" );
	if ( $allow_partial ) {

		if ( $$gene{start_truncated} ) {
			if ( $$gene{orientation} == 1 ) {
				$location = "<$location";
			}
			else {
				$location = ">$location";
			}
		}

		if ( $$gene{stop_truncated} ) {
			my @tmp = split /\.\./, $location;
			my $dnaend = pop @tmp;
			if ( $$gene{orientation} == 1 ) {
				$dnaend = ">$dnaend";
			}
			else {
				$dnaend = "<$dnaend";
			}
			push @tmp, $dnaend;
			$location = join( "..", @tmp );
		}
	}

	return $location;
}

sub get_matpep_location {
	my ( $gene, $pep ) = @_;

	if ( $$gene{is_mutation} ) { return undef }

	my ( $poly_start, $poly_end ) = ( $$pep{pep_start}, $$pep{pep_end} );
	my $location =
		pep_coords_to_dna( $$gene{orientation}, $$gene{exons}, $poly_start, $poly_end, 1 );

	my $allow_partial = get_parameter( "allow_partial_genes" );

	if ( $$pep{start_fuzzy} ) {
		if ( $$gene{orientation} == 1 ) {
			$location = "<$location";
		}
		else {
			$location = ">$location";
		}
	}

	if ( $$pep{end_fuzzy} ) {
		my @tmp = split /\.\./, $location;
		my $dnaend = pop @tmp;
		if ( $$gene{orientation} == 1 ) {
			$dnaend = ">$dnaend";
		}
		else {
			$dnaend = "<$dnaend";
		}
		push @tmp, $dnaend;
		$location = join( "..", @tmp );
	}

	return $location;
}

sub pep_coords_to_dna {
	my ( $ori, $exons, $pep_start, $pep_end, $formatted ) = @_;

	my $cdna_start = 3 * ( $pep_start - 1 ) + 1;
	my $cdna_end = 3 * $pep_end;

	return cdna_coords_to_dna( $ori, $exons, $cdna_start, $cdna_end, $formatted );
}

sub cdna_coords_to_dna {
	my ( $ori, $exons, $cdna_start, $cdna_end, $formatted ) = @_;

	my $dbg = 0;

	my @tmpexons;
	for my $exon ( @{$exons} ) {
		if ( ! exists $$exon{missing_exon} ) {
			push @tmpexons, $exon;
			if ( $dbg ) { print "exon $$exon{dna_start}-$$exon{dna_end}|$$exon{cdna_start}-$$exon{cdna_end}\n" }
		}
	}

	my $exon  = shift @tmpexons;
	while ( $$exon{cdna_end} < $cdna_start && @tmpexons ) {
		if ( $dbg ) { print "skip exon $$exon{dna_start}-$$exon{dna_end}|$$exon{cdna_start}-$$exon{cdna_end}\n" }
		$exon = shift @tmpexons;
	}

	my $dna_start = $$exon{dna_start} + $ori * ( $cdna_start - $$exon{cdna_start} );
	if ( $dbg ) { print " $dna_start = $$exon{dna_start} + $ori * ( $cdna_start - $$exon{cdna_start} )\n" }
	my $formatted_range = "$dna_start..";

	while ( $$exon{cdna_end} < $cdna_end && @tmpexons ) {
		my $exon_end = $$exon{dna_start} + $ori * ( $$exon{cdna_end} - $$exon{cdna_start} );
		$formatted_range .= $exon_end;
		$exon = shift @tmpexons;
		$formatted_range .= ",$$exon{dna_start}..";
		if ( $dbg ) { print "exon $$exon{dna_start}-$$exon{dna_end}|$$exon{cdna_start}-$$exon{cdna_end}\n" }
		if ( $dbg ) { print "$exon_end = $$exon{dna_start} + $ori * ( $$exon{cdna_end} - $$exon{cdna_start} )\n" }
	}
	my $dna_end = $$exon{dna_start} + $ori * ( $cdna_end - $$exon{cdna_start} );
	if ( $dbg ) { print "$dna_end = $$exon{dna_start} + $ori * ( $cdna_end - $$exon{cdna_start} )\n" }
	$formatted_range .= $dna_end;

	if ($formatted) {
		return $formatted_range;
	}
	else {
		return ( $dna_start, $dna_end );
	}
}

sub format_evalue {
	my ( $evalue ) = @_;

	if ( defined $evalue ) {
		my ( $a, $b ) = split /[eE]\-/, $evalue;
		if ( !defined $b || ! length( $b ) ) { $b = 0 }
		if ( $a < 0.0 ) { $a = 0.0 };
		if ( $a > 0 ) {
			while ( $a < 1.0 ) {
				$a = 10. * $a;
				$b++;
			}
			while ( $a >= 10.0 ) {
				$a = $a / 10.0;
				$b--;
			}
			$a = int( 10.0 * $a + 0.5 ) / 10.;
			if ( $a == 10.0 ) {
				$a = 1.0;
				$b++;
			}
			if ( $b >= 1 ) {
				$b =~ s/^0*//;
				$b = lpad( $b, 3, "0" );
			} else {
				$b = "";
			}
		}
		if ( $b ) {
			if ( index( $a, "." ) < 0 ) { $a .= ".0" }
			$evalue = "$a-e$b";
		} else {
			$evalue = $a;
		}
	}

	return $evalue;
}

#
# right pad string to fixed length
sub rpad {
	my ( $text, $pad_len, $pad_char) = @_;

	if ( $pad_len<=0 ) {
		return "";
	} elsif ( $pad_len<=length($text) ) {
		return substr($text,0,$pad_len);
	}

	if ( !defined $pad_char ) {
		$pad_char = " ";
	} elsif ( length($pad_char)>1 ) {
		$pad_char = substr($pad_char,0,1);
	}

	if ( $pad_len>length($text) ) {
		$text .= $pad_char x ( $pad_len - length( $text ) );
	}

	return "$text";
}

#
# left pad string to fixed length
sub lpad {
	my ( $text, $pad_len, $pad_char) = @_;

	if ( $pad_len<=0 ) {
		return "";
	} elsif ( $pad_len<length($text) ) {
		return substr($text,0,$pad_len);
	}

	if ( !defined $pad_char ) {
		$pad_char = " ";
	} elsif ( length($pad_char)>1 ) {
		$pad_char = substr($pad_char,0,1);
	}

	if ( $pad_len>length($text) ) {
		$text = $pad_char x ( $pad_len - length( $text ) ). $text;
	}

	return "$text";
}

#
# center string in fixed length
sub cpad {
	my ( $text, $pad_len, $pad_char) = @_;

	if ( $pad_len<=0 ) {
		return "";
	} elsif ( $pad_len<length($text) ) {
		return substr($text,0,$pad_len);
	}

	if ( !defined $pad_char ) {
		$pad_char = " ";
	} elsif ( length($pad_char)>1 ) {
		$pad_char = substr($pad_char,0,1);
	}

	my $margin = int( ( $pad_len - length($text) ) / 2. );
	if ( $margin>0 ) {
		$text = &lpad($pad_char,$margin,$pad_char) . $text;
	}

	$margin = $pad_len - length($text);
	if ( $margin>0 ) {
		$text .= &rpad($pad_char,$margin,$pad_char);
	}

	return "$text";
}

sub print_hash {
	my ( $hashname, $hash ) = @_;

	our %hashes;
	#if ( exists $hashes{$hashname} ) { return }
	$hashes{$hashname} = $hash;

	print "\n$hashname\n";
	for my $hashkey ( sort keys %$hash ) {
		if ( !defined $$hash{$hashkey} ) {
			print "$hashkey=\n";
			next;
		}
		my $hashval = "$$hash{$hashkey}";
		if ( $hashval =~ /HASH/ ) {
			print_hash( "$hashname.$hashkey", $$hash{$hashkey} );
		}
		elsif ( $hashval =~ /ARRAY/ ) {
			for my $element ( @{ $$hash{$hashkey} } ) {
				if ( "$element" =~ /HASH/ ) {
					print_hash( "$hashname.$hashkey.element", $element );
				}
				else {
					print "$hashname=>ARRAY: $$hash{$hashkey}\n";
				}
				last;
			}
		}
		else {
			print "$hashkey=$hashval\n";
		}
	}
	return;
}

#==========================================================================
# customization stubs

sub customize_defaults {
	return;
}

# ---- FIND_CANDIDATES ----
# customized initialization of find_candidates
# (e.g. tweak parameters for specific genomic sequence)
sub before_find_candidates {
	my ( $genome ) = @_;
	return;
}

# override default comparison of candidate regions
# return1 to prefer ihit, -1 to prefer jhit, 0 to use default comparison
sub candidate_comparison_exception {
	my ( $genome, $ihit, $jhit ) = @_;
	return undef;
}

# customized post processing of candidate regions
sub after_find_candidates {
	my ( $genome, $candidates ) = @_;
}


# ---- FIND_RAWGENES ----
# customized initialization of find_rawgenes
# (e.g. tweak parameters for specific genomic sequence)
sub before_find_rawgenes {
	my ( $genome, $candidate ) = @_;
	return;
}

# override default selection (highest bit score) of raw genes
# return array reference
sub select_rawgenes_exception {
	my ( $genome, $candidate, $rawgenes ) = @_;
	return undef;
}

# customized post-processing of raw genes
sub after_find_rawgenes {
	my ( $genome, $candidate, $rawgenes ) = @_;
	return;
}

# when multiple permutations exist for a rawgene
# allow the user to select among them
sub screen_rawgene_permutations {
	my ( $genome, $rawgene ) = @_;
	return;
}

# ---- FIND_GENE ----
# customized initialization of find_gene
# (e.g. tweak parameters for specific genomic sequence)
sub before_find_gene {
	my ( $genome, $rawgene ) = @_;
	return;
}

# override default selection of stop codon
# return genomic position of last base of stop codon (or of last peptide if partial)
sub stop_codon_exception {
	my ( $sequence, $rawgene, $exons ) = @_;
	return undef;
}

# override default selection of start codon
# return genomic position of first base of start codon (or of first peptide if partial)
sub start_codon_exception {
	my ( $sequence, $rawgene, $exons ) = @_;
	return undef;
}

# customized post-processing of gene
sub after_find_gene {
	my ( $genome, $rawgene, $gene ) = @_;
	return;
}

# customized validation of gene
sub before_validate_gene {
	my ( $genome, $tmpgene ) = @_;
	return undef;
}

sub after_validate_gene {
	my ( $genome, $tmpgene ) = @_;
	return undef;
}

# split polyprotein into mature peptides
# return array reference
sub get_mature_peptides {
	my ($gene) = @_;
	return undef;
}

# customize the genomic sequence's set of genes
sub customize_gene_set {
	my ( $genome, $genes ) = @_;
	return;
}

# --- REPORT_GENE ----
# add stuff before standard report
sub before_report_gene {
	my ( $TBL, $CDS, $PEP, $ALIGN, $FS, $genome, $genes ) = @_;
	return;
}

# exceptions to standard report
sub report_gene_exception {
	my ( $TBL, $CDS, $PEP, $ALIGN, $FS, $genome, $genes ) = @_;
	return ( 1, 1, 1, 1 );
}

# add stuff after standard report
sub after_report_gene {
	my ( $gene, $TBL, $CDS, $PEP, $ALIGN, $FS ) = @_;
	return;
}

# ---- standardize reference data ----
# return standardized reference name (name is NOT required)
sub get_reference_name {
	my ($reference_id) = @_;

	my $seq = get_reference_seq($reference_id);
	if ( !defined $seq ) { return $reference_id }

	if ( " $$seq{defline} " =~ / gene="([^"]+)"/i ) {
		my $name = $1;
		$name =~ s/^ +//;
		$name =~ s/ +$//;
		return $name;
	}
	elsif ( " $$seq{defline} " =~ / gene=([^ ]+) /i ) {
		my $name = $1;
		$name =~ s/^ +//;
		$name =~ s/ +$//;
		return $name;
	}
	elsif ( $$seq{defline} =~ /^ *([^ ].{1,8}) *$/ ) {
		my $name = $1;
		return $name;
	}
	elsif ( $$seq{defline} =~ /poly *protein/i ) {
		my $name = "POL";
		if ( "$$seq{defline} " =~ /(gene|orf) *([0-9])[ -]*([ab]*)[\W]/i ){
			$name .= $1.$2;
		}
		return $name;
	}
	else {
		return $reference_id;
	}
}

# return standardized reference definition (definition IS required)
sub get_reference_definition {
	my ($reference_id) = @_;

	my $seq = get_reference_seq($reference_id);
	if ( !defined $seq ) { return undef }

	my $definition = $$seq{defline};
	$definition =~ s/\s+/ /g;
	$definition =~ s/^ +//;
	$definition =~ s/^[^ ]* *//;
	$definition =~ s/ +$//;

	if ( " $definition " =~ / product="([^"]+)"/i ) {
		$definition = $1;
		$definition =~ s/^ +//;
		$definition =~ s/ +$//;
	}
	elsif ( " $definition " =~ / product=([^ ]+) /i ) {
		$definition = $1;
		$definition =~ s/^ +//;
		$definition =~ s/ +$//;
	}

	return $definition;
}

sub check_coverage {
	my ( $min_coverage, $feature, $parent, $dbg ) = @_;
	if ( ! defined $dbg ) { $dbg = 0 }

	my $coverage;
	if ( exists $$feature{subject_id} ) {
		$coverage = $$feature{pct_scoverage};
	}
	else {
		$coverage = $$feature{pct_refcoverage};
	}
	if ( $coverage >= $min_coverage ) {
		if ( $dbg ) { print "$coverage >= $min_coverage\n" }
		return 1;
	}
	if ( ! get_parameter( "allow_partial_genes" ) ) {
		if ( $dbg ) { print "$coverage < $min_coverage | no partial genes\n" }
		return 0;
	}

	my $qleft;
	my $qright,
	my $ori;
	my $origcoverage = $coverage;
	my $missing5;
	my $missing3;
	my $spliced;

	if ( exists $$feature{subject_id} ) {
		$qleft = $$feature{query_left};
		$qright = $$feature{query_right};
		$ori = $$feature{orientation};
		$missing5 = 100.0 * ( $$feature{subject_left} - 1 ) / $$feature{subject_length};
		$missing3 = 100.0 * ( $$feature{subject_length} - $$feature{subject_right} ) / $$feature{subject_length};
	}
	else {
		$ori = $$feature{orientation};
		if ( $ori == 1 ) {
			$qleft = $$feature{start_codon};
			$qright = $$feature{stop_site};
		}
		else {
			$qleft = $$feature{stop_site};
			$qright = $$feature{start_codon};
		}
		$missing5 = $$feature{pct_reftrunc5};
		$missing3 = $$feature{pct_reftrunc3};
	}
	if ( $coverage + $missing5 + $missing3 < $min_coverage ) {
		if ( $dbg ) { print "$coverage + $missing5 + $missing3 < $min_coverage\n" }
		return 0;
	}

	my $max_distance;
	my $parent_length;
	if ( exists $$parent{sequence} ) {
		$parent_length = length( $$parent{sequence} );
		$max_distance = 45;
		if ( allow_splicing($feature) ) { $max_distance += get_parameter( "max_intron_size" ) + 20 }
	}
	else {
		$parent_length = $$parent{protein_length};
		$max_distance = 15;
	}

	my $maxleft = $max_distance;
	my $minright = $parent_length - $max_distance + 1;

	if ( $missing5 > 0 ) {
		if ( $ori == 1 ) {
			if ( $qleft <= $maxleft ) { $coverage += $missing5 }
		}
		else {
			if ( $qright >= $minright ) { $coverage += $missing5 }
		}
	}

	if ( $missing3 > 0 ) {
		if ( $ori == 1 ) {
			if ( $qright >= $minright ) { $coverage += $missing3 }
		}
		else {
			if ( $qleft <= $maxleft ) { $coverage += $missing3 }
		}
	}

	if ( $dbg ) { print "original=$origcoverage  missing5/3=$missing5/$missing3  range=$maxleft-$minright  feature=$qleft-$qright  modified coverage=$coverage\n" }

	if ( $coverage >= $min_coverage ) { return 1 }
	return 0;
}

sub parse_exception_lists {

	my $list = get_parameter( "spliced_gene_list" );
	for my $gene ( split /,/, $list ) {
		$splicedGenes{$gene} = 1;
	}

	$list = get_parameter( "ribosomal_slippage_list" );
	for my $gene ( split /,/, $list ) {
		$slippageGenes{$gene} = 1;
	}

	$list = get_parameter( "stop_readthru_list" );
	for my $gene ( split /,/, $list ) {
		$readthruGenes{$gene} = 1;
	}
}

sub allow_stopcodon_readthru {
	my ( $hit ) = @_;

	my $allow = get_parameter( "allow_stopcodon_readthru" );
	if ( $allow == 0 || $allow == 1) { return $allow }

	my $defline = $$hit{subject_definition};
	if ( ! defined $defline ) {
		$defline = "";
		if ( exists $$hit{ref_id} ) {
			my $ref = get_reference_seq( $$hit{ref_id} );
			$defline = $$ref{defline};
		}
	}
	if ( " $defline " =~ / stopcodon_readthru="{0,1}Y/ ) {
		return 1;
	}


	my $name;
	if ( exists $$hit{subject_id} ) {
		$name = get_reference_name($$hit{subject_id});
	}
	else {
		$name = get_reference_name($$hit{ref_id})
	}
	if ( ! defined $name ) {
		return 0;
	}

	if ( exists $readthruGenes{$name} ) {
		return 1;
	}
	else {
		return 0;
	}
}

sub allow_ribosomal_slippage {
	my ( $hit ) = @_;

	my $allow = get_parameter( "allow_ribosomal_slippage" );
	if ( $allow == 0 || $allow == 1) { return $allow }

	my $defline = $$hit{subject_definition};
	if ( ! defined $defline ) {
		$defline = "";
		if ( exists $$hit{ref_id} ) {
			my $ref = get_reference_seq( $$hit{ref_id} );
			$defline = $$ref{defline};
		}
	}
	if ( " $defline " =~ / ribosomal_slippage="{0,1}Y/ ) {
		return 1;
	}

	my $name;
	if ( exists $$hit{subject_id} ) {
		$name = get_reference_name($$hit{subject_id});
	}
	else {
		$name = get_reference_name($$hit{ref_id})
	}
	if ( ! defined $name ) {
		return 0;
	}

	if ( exists $slippageGenes{$name} ) {
		return 1;
	}
	else {
		return 0;
	}
}

sub allow_splicing {
	my ( $hit ) = @_;
	my $allow = get_parameter( "allow_splicing" );
	if ( $allow == 0 || $allow == 1) { return $allow }

	my $defline = $$hit{subject_definition};
	if ( ! defined $defline ) {
		$defline = "";
		if ( exists $$hit{ref_id} ) {
			my $ref = get_reference_seq( $$hit{ref_id} );
			$defline = $$ref{defline};
		}
	}
	if (  " $defline " =~ / spliced="{0,1}Y/i ) {
		return 1;
	}

	my $name;
	if ( exists $$hit{subject_id} ) {
		$name = get_reference_name($$hit{subject_id});
	}
	else {
		$name = get_reference_name($$hit{ref_id})
	}
	if ( ! defined $name ) {
		return 0;
	}

	if ( exists $splicedGenes{$name} ) {
		return 1;
	}
	else {
		return 0;
	}
}
1;
