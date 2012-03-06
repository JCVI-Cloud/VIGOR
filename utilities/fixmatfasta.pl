#!/usr/local/bin/perl
use strict;

my @fasta = split /\n/, `cat Rhinovirus_mature_pep_db.fasta.new`;

my %defs;
for my $line ( split /\n/, `cat matdat.tsv` ) {
	my ( $id, $def ) = split /\t/, $line;
	$defs{$id} = $def;
print "$id\t$def\n";
}

open( FASTA, ">Rhinovirus_mature_pep_db.fasta" );
for my $line ( @fasta ) {
	if ( $line =~ /^>/ ) {
		$line =~ s/ .*$//;
		$line .= " " . $defs{substr($line,1)};		
	}
	print FASTA "$line\n";
	print "$line\n";
}
close FASTA;
