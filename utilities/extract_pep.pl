#!/usr/local/bin/perl
use strict;
use Getopt::Std;
use Cwd 'realpath';
use File::Basename;

my $program = realpath($0);
my $myBin = dirname($program);
require "$myBin/VIGOR.pm";

our ( $opt_i );
getopt('i');

my %peps;
print "loading $opt_i\n";
my %seqs = loadFasta( $opt_i ); 
for my $seq ( sort { $$a{id} cmp $$b{id} } values %seqs ) {
	print "$$seq{defline}\n";
	if ( $$seq{defline} =~ / definition=\".*(nsp[0-9]+).*\"/ ) {
		my $nsp = lc $1;
		$$seq{defline} =~ s/\.[^ ]* /_$nsp /;
		$$seq{id} =~ s/\..*$/_$nsp/;
		$peps{$nsp}{$$seq{id}} = $seq;
		print "** $nsp: $$seq{defline}\n";
	}
	elsif ( $$seq{defline} =~ / gene=([^ ]+) / ) {
		my $gene = $1;
		if ( $gene eq "NSP" ) {
			$gene = "3";
			$$seq{defline} =~ s/ gene=NSP / gene=3 /;
		}
		$$seq{defline} =~ s/\.[^ ]* /$gene /;
		$$seq{id} =~ s/\..*$/$gene/;
		$peps{$gene}{$$seq{id}} = $seq;
	}
}

for my $pep ( sort { substr( $a, 3 ) <=> substr( $b, 3 ) } keys %peps ) {
	open ( PEP, ">$opt_i" . "_$pep" );
	for my $id ( sort keys %{$peps{$pep}} ) {
		my $seq = $peps{$pep}{$id};
		print PEP ">$$seq{defline}\n";
		$$seq{sequence} =~ s/(.{60})/$1\n/g;
		if ( substr( $$seq{sequence}, length( $$seq{sequence} ) - 1, 1 ) ne "\n" ) {
			$$seq{sequence} .= "\n";
		}
		print PEP $$seq{sequence};
	}
	close PEP;
}

exit(0);