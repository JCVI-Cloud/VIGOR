#!/usr/local/bin/perl -w
use strict;
use Getopt::Std;
use Cwd 'realpath';
use File::Basename;

my $program = realpath($0);
our $myBin = dirname($program);
our $myData = "$myBin/data";
require "$myBin/VIGOR.pm";

our ( $opt_g );
getopt('g');

*STDERR = *STDOUT;
$|++;


my %genes;
my %proteins;
my $acc;
my $pepno;
my $pepid;
my $geneid;
my $product = "predicted protein";
my $spliced = "N";
my $sequence;
while ( my $line = <STDIN> ) {
	chomp $line;
	$line =~ s/\t/        /g;

#print "$line\n";
	if ( $line =~ /^LOCUS/ ) {
		( undef, $acc ) = split / +/, $line;
		$pepno = 0;
	}
	elsif ( $line =~ /^ * (CDS|gene) .*join/i ) {
		$spliced="Y";
	}
	elsif ( $line =~ /^ *\/gene="([^"]+)"/ ) {
		$geneid = $1;
	}
	elsif ( $line =~ /^ *\/protein_id="([^"]+)"/ ) {
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
		$pepno++;
		$sequence =~ s/[^a-zA-Z]//g;
#print "sequence=$sequence\n";
		$sequence =~ s/(.{60})/$1\n/g;
		if ( ! defined $geneid ) {
			if ( $opt_g ) {
				$spliced = "N";
				$geneid = undef;
				$pepid = undef;
				$product = "predicted protein";
				$sequence = undef;
				next;
			}
			$geneid = "$acc-g$pepno";
		}
		if ( ! defined $pepid ) {
			$pepid = "$acc-p$pepno";
		}
		if ( exists $proteins{$pepid} ) {
			my $subid = 2;
			while ( exists $genes{"$pepid-p$subid"} ) {
				$subid++;
			}
			$pepid = "$pepid-p$subid";
		}
		if ( exists $genes{$geneid} ) {
			my $subid = 2;
			while ( exists $genes{"$geneid-p$subid"} ) {
				$subid++;
			}
			$geneid = "$geneid-p$subid";
		}
		$genes{$geneid} = $pepid;
		$proteins{$pepid} = $geneid;

		print ">$pepid gene=\"$geneid\" product=\"$product\" spliced=$spliced\n$sequence\n";
		$spliced = "N";
		$geneid = undef;
		$pepid = undef;
		$product = "predicted protein";
		$sequence = undef;
	}
}
close STDIN;
exit(0);
