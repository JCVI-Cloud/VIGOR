#!/usr/local/bin/perl
use strict;
use Getopt::Std;
use Cwd 'realpath';
use File::Basename;

my $program = realpath($0);
my $myBin = dirname($program);
require "$myBin/VIGOR.pm";

our ( $opt_t );
getopt('t');
my $taxonomy = lc $opt_t;
if ( !defined $taxonomy ) { $taxonomy = "viruses;" }

my $fileprefix = $taxonomy;
$fileprefix =~ s/[^a-z0-9A-Z_]//g;
open( PEP, ">$fileprefix.mature_pep.faa" );
open( CDNA, ">$fileprefix.mature_pep.fna" );

my @genomes = split /\n/, `grep -i "$taxonomy" /usr/local/projects/CAMERA/runtime-shared/filestore/system/GenomeProject/allGenomes.tsv`;
for my $genome ( @genomes ) {
	my ( @tmp ) = split /\t/, $genome;
	my $dir = pop @tmp;
	my @files = split /\n/, `ls -1 $dir/*.gbk`;
	for my $file ( @files ) {
print "$file\n";		
		my $acc = `grep VERSION $file | head -n 1`;
		( $acc ) = split /\n/, $acc;
		$acc =~ s/^ *VERSION *//;
		$acc =~ s/ .*$//;
#print "ACC=$acc\n";
		
		my $org = `grep ORGANISM $file | head -n 1`;
		( $org ) = split /\n/, $org;
		$org =~ s/^ *ORGANISM *//;
#print "ORG=$org\n";

		my $gbk = `cat $file`;
		my @pieces = split /ORIGIN/, $gbk;
		my @tmpseq = split /\n/, pop @pieces;
		my $seq = join( "", @tmpseq );
		$seq =~ s/[^a-zA-Z]//g;
#print "SEQ=$seq\n";

		$gbk = join( "\t", split /\n/, join( "ORIGIN", @pieces ) );
		@pieces = split /(  mat_peptide  )/i, $gbk;
		if ( @pieces < 2 ) { next }

		my $mpno = 0;
		my $p = 0;
		while ( $p < @pieces ) {
			my $piece = $pieces[$p];
			if ( $piece =~ /  mat_peptide  /i ) {
#print "$piece\n";
				$mpno++;
				
				$p++;
				$piece = $pieces[$p];
				my @mptmp = split /\t/, $piece;
				my $position = shift @mptmp;
				$position =~ s/ //g;
#print "  $position\n";
				my $coords = $position;
				$coords =~ s/[[\]{}<>()]//g;
				$coords =~ s/join//gi;

				my $data ="";
				$p++;
				for my $tmp ( @mptmp ) {
					if ( $tmp !~ /^      / ) { last }
#print "  $tmp\n";
					$tmp =~ s/^ *//;
					$data .= " " . $tmp;
				}
				
				my $ori = 1;
				if ( $coords =~ /reverse/i || $coords =~ /complement/i ) {
					$ori = -1;
					$coords =~ s/reverse//gi;
					$coords =~ s/complement//gi;
				}
				
				my @segments;
				for my $coord ( split /,/, $coords ) {
					my ( $left, $right ) = split /\.\./, $coord;
					my %segment;
					if ( $ori == 1 ) {
						$segment{start} = $left;
						$segment{end} = $right;
					}
					else {
						$segment{start} = $right;
						$segment{end} = $left;
					}
					push @segments, \%segment;
				}
				if ( $ori == -1 ) {
					@segments = reverse @segments;
				}
				
				my $cdna = "";
				for my $segment( @segments ) {
					my $exon = subsequence( $seq, $$segment{start}, $$segment{end} );
					$cdna .= $exon;
				}
				my $pep = DNA2AA( $cdna );
				my $cdnalen = length( $cdna );
				my $peplen = length( $pep );
				
				my $protid = $acc . "_mp" . $mpno;
				if ( $data =~ /\/protein_id=\"(.*)\"/ ) {
					$protid = $1;
					$data =~ s/ \/protein_id=\".*\"*//;
				}

				my $def = ">$protid /genome_acc=\"$acc\" /organism=\"$org\" /position=\"$position\" /protein_length=$peplen" . $data;
				print "  $def\n";
				
				$pep =~ s/(.{60})/$1\n/g;
				if ( substr( $pep, length( $pep ) - 1, 1 ) ne "\n" ) { $pep .= "\n" } 
				print PEP "$def\n$pep";
				
				$cdna =~ s/(.{60})/$1\n/g;
				if ( substr( $cdna, length( $cdna ) - 1, 1 ) ne "\n" ) { $cdna .= "\n" } 
				print CDNA "$def\n$cdna";
			}
			else {
				$p++;				
			}
		}
		
	}
}
close PEP;
close CDNA;		
 