#!/usr/local/bin/perl -w
use strict;
use Getopt::Std;
use Cwd 'realpath';
use File::Basename;

my $program = realpath($0);
my $myBin = dirname($program);
require "$myBin/VIGOR.pm";

my %na = loadFasta( $ARGV[0] );
open( OUT, ">$ARGV[1]" );

for my $id ( sort keys %na ) {
	my $aa = $na{$id};
	$$aa{sequence} = DNA2AA( $$aa{sequence} );
	$$aa{sequence} =~ s/(.{60})/$1\n/g;
	if ( substr( $$aa{sequence}, length($$aa{sequence})-1, 1 ) ne "\n" ) { $$aa{sequence} .= "\n" }
	print OUT ">$$aa{defline}\n$$aa{sequence}";
}

close OUT;
exit(0);
