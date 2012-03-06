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

my %nsps;
print "loading $opt_i\n";
my %seqs = loadFasta( $opt_i ); 
for my $seq ( sort { $$a{id} cmp $$b{id} } values %seqs ) {
	print "$$seq{defline}\n";
	if ( $$seq{defline} =~ / definition=\".*(nsp[0-9]+).*\"/ ) {
		my $nsp = lc $1;
		$$seq{defline} =~ s/\.[^ ]* /_$nsp /;
		$$seq{id} =~ s/\..*$/_$nsp/;
		$nsps{$nsp}{$$seq{id}} = $seq;
		print "** $nsp: $$seq{defline}\n";
	}
}

for my $nsp ( sort { substr( $a, 3 ) <=> substr( $b, 3 ) } keys %nsps ) {
	open ( NSP, ">$opt_i" . "_$nsp" );
	for my $id ( sort keys %{$nsps{$nsp}} ) {
		my $seq = $nsps{$nsp}{$id};
		print NSP ">$$seq{defline}\n";
		$$seq{sequence} =~ s/(.{60})/$1\n/g;
		if ( substr( $$seq{sequence}, length( $$seq{sequence} ) - 1, 1 ) ne "\n" ) {
			$$seq{sequence} .= "\n";
		}
		print NSP $$seq{sequence};
	}
	close NSP;
}

exit(0);