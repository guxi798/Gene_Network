#!/usr/bin/perl
use strict;
use warnings;

my $srcfolder = $ARGV[0];
my $tgtfile = $ARGV[1];

opendir(SRC, "$srcfolder") or die "Cannot open $srcfolder: $!";
open(TGT, ">$srcfolder/$tgtfile") or die "Cannot write $tgtfile: $!";

my @subs = sort(grep(/Gsea/, readdir(SRC)));
foreach my $sub (@subs) {
	opendir(SUB, "$srcfolder/$sub") or die "Cannot open $sub: $!";
	my @files = sort {$b cmp $a} grep(/^gsea_report.+xls$/, readdir(SUB));
	foreach my $file (@files) {
		open(FILE, "$file") or die "Cannot open $file: $!";
		foreach my $line (<FILE>) {
			print TGT $line;
		}
		close(FILE);
	}
	close(SUB);
}

close(SRC);
close(TGT);