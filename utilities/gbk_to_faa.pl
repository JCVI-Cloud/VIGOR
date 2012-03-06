#!/usr/local/bin/perl -w
use strict;
use Getopt::Std;
use Cwd 'realpath';
use File::Basename;

my $program = realpath($0);
my $myBin = dirname($program);
require "$myBin/VIGOR.pm";

*STDERR = *STDOUT;
$|++;

my $pepid;
my $product = "predicted protein";;
my $sequence;
while ( my $line = <STDIN> ) {
	chomp $line;
	$line =~ s/\t/        /g;

#print "$line\n";
	if ( $line =~ /^ *\/protein_id="([^"]+)"/ ) {
		$pepid = $1;
#print "pepid=$pepid\n";
	}
	elsif ( $line =~ /^ *\/product="([^"]+)"/ ) {
		$product = $1;
	}
	elsif ( $line =~ /^ *\/translation="([^"]+)/ ) {
		$sequence = $1;
		while ( my $line = <STDIN> ) {
			chomp $line;
#print "$line\n";
			$sequence .= $line;
			if ( $line =~ /" *$/ ) { last }
		}
		$sequence =~ s/[^a-zA-Z]//g;
#print "sequence=$sequence\n";
		$sequence =~ s/(.{60})/$1\n/g;
	}
	elsif ( $line =~ /^ *gene  /i || $line =~ /^ *CDS   / ) {
		if ( defined $sequence ) { print ">$pepid $product\n$sequence\n" }	
		$pepid = undef;
		$product = "predicted protein";
		$sequence = undef;
	}
}
if ( defined $pepid ) { print ">$pepid $product\n$sequence\n" }		
close STDIN;
exit(0);
