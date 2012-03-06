#!/usr/local/bin/perl -w
use strict;
use Getopt::Std;
use Cwd 'realpath';
use File::Basename;

my $program = realpath($0);
my $myBin = dirname($program);
require "$myBin/VIGOR.pm";

my %vdb = loadFasta("/usr/local/projects/PANDA/db/panda/ViralGroup/ViralGroup.niaa");
for my $seq ( sort { $$a{id} cmp $$b{id} } values %vdb ) {
	if ( $$seq{defline} =~ /1 *ab/i ) { next };
	if ( $$seq{defline} =~ /^PDB/ ) { next };
	if ( $$seq{defline} =~ /^PIR/ ) { next };
	if ( $$seq{defline} =~ /^OMNI/ ) { next };
	if ( $$seq{defline} =~ / *taxon:[^{]*\{/ ) {
		$$seq{defline} =~ s/ *taxon:[^{]*(\{[^{]*\}).*$/ $1/;
	}
	else {
		$$seq{defline} =~ s/ *taxon:.*$//;
	}
	$$seq{defline} =~ s/;\}/}/;
	$$seq{defline} =~ s/^[^ ]*  *//;
	
	$$seq{sequence} =~ s/(.{60})/$1\n/g;
	if ( substr( $$seq{sequence}, length( $$seq{sequence} ) - 1, 1 ) ne "\n" ) {
		$$seq{sequence} .= "\n";
	}
	print ">$$seq{id} $$seq{defline}\n$$seq{sequence}";
}
exit(0);