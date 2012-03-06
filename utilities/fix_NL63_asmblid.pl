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

print "loading asmbl_id2sample_id.txt\n";
my $fixfile = `cat asmbl_id2sample_id.txt`;
my %fixes;
for my $fixline ( split /\n/, $fixfile ) {
	my ( $asmblid, $sampleid ) = split /\t/, $fixline;
	$fixes{$asmblid} = $sampleid;
}

my %nsps;
print "loading $opt_i\n";
my %seqs = loadFasta( $opt_i ); 
open( OUT, ">$opt_i.fixed" );
for my $seq ( sort{ $$a{id} <=> $$b{id} } values %seqs ) {
	my $id = $$seq{id};
	if ( exists $fixes{$$seq{id}} ) {
		$id = "$fixes{$$seq{id}}:$id";
	}
print "\"$$seq{id}\" => $id\n";
	my $defline = $$seq{defline};
	$defline =~ s/^[^ ]* *//;
	$defline = ">$id $defline";
	my $seq = $$seq{sequence};
	$seq =~ s/(.{60})/$1\n/g;
	if ( substr( $seq, length($seq)-1, 1 ) ne "\n" ) { $seq .= "\n" } 
	print OUT "$defline\n$seq";
}
close OUT;
exit(0);