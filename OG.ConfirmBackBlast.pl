###########################################################
###	@Author		Xi Gu									###
###	@Time		July 9, 2014							###
###	@Function	Confirm OG Results						###
###########################################################

#!/usr/bin/perl
use strict;

###########################################################################################
###  Input arguments. As described in the Input section above.							###
###########################################################################################

my $sourcefile = "../14.backblast/blast.MBH.tab";								## path to the directory containing all proteomes
my $sourcefile2 = "../13.orthomclMclToGroups/phytozome.singletons.txt";
my $outfile = "../14.backblast/OG.MBH.tab";									## path to the target directory
my $ogfile = "../13.orthomclMclToGroups/ortholog_groups.tab";		## ortholog group file produced by orthoMCL

###########################################################################################
###########################################################################################

my %hash = ();													## initialize the big hash

open(SRC, "$sourcefile") or die "Cannot open $sourcefile: $!";

foreach my $line (<SRC>) {
	chomp $line;
	my @lines = split(/\t/, $line);

	if(not exists $hash{$lines[0]}){
		$hash{$lines[0]} = [$lines[1]];
	}
	else{
		push @{$hash{$lines[0]}}, $lines[1];
	}
}

close(SRC);

my %singleton = ();
open(SRC, "$sourcefile2") or die "Cannot open $sourcefile2: $!";

foreach my $line (<SRC>) {
	chomp $line;
	my @lines = split(/[|\.]/, $line);
	my $sp = shift(@lines);
	pop(@lines);
	my $gene = join(".", @lines);

	if(not exists $singleton{$gene}){
		$singleton{$gene} = 0;
	}

}

close(SRC);


###########################################################################################
###########################################################################################

open(SRC, "$ogfile") or die "Cannot open $ogfile: $!";				## read in ortholog group file
open(TGT, ">$outfile") or die "Cannot write $outfile: $!";

foreach my $line (<SRC>) {
	chomp $line;
	my @lines = split(/\t/, $line);							## get all genes in the target OG
	my $group = shift(@lines);								## group name
	my $miss_sp = 0;
	my $count_in = 0;
	my $count_out = 0;

	if ($line !~ /AT/ and $line =~ /Glyma/ and $line =~ /Potri/ and $line =~ /Phvul/){
		$miss_sp = "ath";
	}elsif($line =~ /AT/ and $line !~ /Glyma/ and $line =~ /Potri/ and $line =~ /Phvul/){
		$miss_sp = "gma";
	}elsif($line =~ /AT/ and $line =~ /Glyma/ and $line !~ /Potri/ and $line =~ /Phvul/){
		$miss_sp = "ptr";
	}elsif($line =~ /AT/ and $line =~ /Glyma/ and $line =~ /Potri/ and $line !~ /Phvul/){
		$miss_sp = "pvu";
	}else{
		next;
	}

	print TGT "$group\t";
	my $single = "FALSE";
	my %temp = ();
	foreach my $gene (@lines) {
		$count_out++;
		if(exists $hash{$gene}){
			my $prefix = &species2prefix($miss_sp);
			my $mark = 0;
			foreach my $ortho (@{$hash{$gene}}) {
				if($ortho =~ /$prefix/){
					$temp{$gene} = $ortho;
					$mark = 1;
					$count_in++;
					$single = "TRUE" if (exists $singleton{$ortho});
				} # if
			} # foreach
			if($mark == 0){
				$temp{$gene} = 0;
			}
		} #if
		else{
			$temp{$gene} = "NA";
		}
	} # foreach
	print TGT "$count_out\t$count_in\t$single\t";
	foreach my $key (sort(keys %temp)) {
		print TGT "$key:",$temp{$key},"\t";
	}
	print TGT "\n";
}
close(TGT);
close(SRC);

###########################################################################################
###  Self-define function.																###
###########################################################################################

sub species2prefix{													## sp name to gene prefix
	my $sp = shift;
	return "AT" if $sp eq "ath";
	return "Glyma" if $sp eq "gma";
	return "Potri" if $sp eq "ptr";
	return "Phvul" if $sp eq "pvu";
}

sub prefix2species{													## gene prefix to sp name
	my $pf = shift;
	return "ath" if $pf =~ "^AT";
	return "gma" if $pf eq "Glyma";
	return "ptr" if $pf eq "Potri";
	return "pvu" if $pf eq "Phvul";
}

