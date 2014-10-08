###########################################################
###	@Author		Xi Gu									###
###	@Time		Aug 05, 2014							###
###	@Function	Assess Annotations						###
###########################################################

###########################################################################################
###  This script assesses the annotation for different species by compare them for each	###
###  ortholog group defined by orthoMCL. It is focused on KEGG and GO databases which	###
###  will be used for downstream GSEA and GO enrichment analysis.						###
###  For KEGG, this annotation mainly comes from Arabidopsis, so we assigned orthologs	###
###  in other species which don't have KEGG terms while their inter-species orthologs	###
###  have. Basically, the operation is to transfer available KEGG terms to unannotated	###
###  orthologs.																			###
###  For GO, we first check to what extent orthologs in the same group but from diff    ###
###  species agree with each other (consistency level) and print the Jaccord coef to	###
###  output files. The output files will be further manually checked but it's out of	###
###  this script's scope.																###
###																						###
###  Input:  1) annotation files of all species, downloaded from phytozome.				###
###			 2) ortholog group file.													###
###	 Output: 1) updated annotation files for all sp with edited KEGG annotation.		###
###			 2) soft links for latest updated annotation files.							###
###			 3) assessment information for KEGG database.								###
###			 4) assessment information for GO database.									###
###########################################################################################

###########################################################################################
###  Load necessary packages.															###
###########################################################################################

#!/usr/bin/perl
use strict;
use Array::Utils qw(:all);								## provide interaction and union
use List::Compare;

###########################################################################################
###  Provide location of annotation files and ortholog group file.						###
###########################################################################################

#my $dir = "/lustre1/escratch1/guxi798_Jul_17";			## home dir
my @sourcefiles = qw(./arabidopsis/00.data/annotation_info_latest.tab
				   ./populus/00.data/annotation_info_latest.tab
				   ./soybean/00.data/annotation_info_latest.tab
				   ./common.bean/00.data/annotation_info_latest.tab);
my @tgtannot = qw(./arabidopsis/00.data/Athaliana_167_TAIR10.annotation_edit_V3.0.txt
				   ./populus/00.data/Ptrichocarpa_210_v3.0.annotation_edit_V3.0.txt
				   ./soybean/00.data/Gmax_189_annotation_edit_V3.0.txt
				   ./common.bean/00.data/Pvulgaris_218_annotation_edit_V3.0.txt);
my @sp = ("ath", "ptr", "gma", "pvu");					## species name

###########################################################################################
###  Use hash structure to store all annotation information. The hash contain 3 levels:	###
###  At the first level, the key is species name; at the second level, the key is gene 	###
###  name; at the third level, the key is annotation database name. The value at leaf	###
###  node is specific annotation.														###
###  Note: this script cannot handle AS transcript annotation.							###
###########################################################################################

print "******************** Construct hash structure to store annotation *******************\n";
my %bighash = ();										## store the whole annotation
foreach my $file (@sourcefiles) {						## for each species
	my %sphash = ();									## sub-hash for each species
	open(SRC, "$file") or die "Cannot open $file: $!";

	foreach my $line (<SRC>) {							## for each gene entry
		chomp $line;
		my @lines = split(/\t/, $line);
		if (not exists $sphash{$lines[1]}) {			## key is gene name
			$sphash{$lines[1]}{"ID"} = [split(/,/, $lines[0])];
			$sphash{$lines[1]}{"TName"} = [split(/,/, $lines[2])];
			$sphash{$lines[1]}{"PName"} = [split(/,/, $lines[3])];
			$sphash{$lines[1]}{"PFAM"} = [split(/,/, $lines[4])];		## database name
			$sphash{$lines[1]}{"Panther"} = [split(/,/, $lines[5])];
			$sphash{$lines[1]}{"KOG"} = [split(/,/, $lines[6])];
			$sphash{$lines[1]}{"EC"} = [split(/,/, $lines[7])];
			$sphash{$lines[1]}{"KEGG"} = [split(/,/, $lines[8])];
			$sphash{$lines[1]}{"GO"} = [split(/,/, $lines[9])];
			$sphash{$lines[1]}{"HitName"} = [split(/,/, $lines[10])];
			$sphash{$lines[1]}{"Symbol"} = [split(/,/, $lines[11])];
			$sphash{$lines[1]}{"HitDefline"} = $lines[12];
		}
	}
	close(SRC);

	my $sp = shift(@sp);								## species name
	$bighash{$sp} = {%sphash};							## key is sp name
}

print "******************** Done *******************\n\n";

###########################################################################################
###  Read in GO assessment result file.	###
###########################################################################################

system("cp ./common.bean/03.result/09.annotAssessment/02.annotAssessment.GO.tab ./common.bean/03.result/09.annotAssessment/02.pre.annotAssessment.GO.tab");
system("cp ./common.bean/03.result/09.annotAssessment/03.annotAssessment.GO.detail.tab ./common.bean/03.result/09.annotAssessment/03.pre.annotAssessment.GO.detail.tab");

my $assfile = "./common.bean/03.result/09.annotAssessment/03.pre.annotAssessment.GO.detail.tab";

open(SRC, "$assfile") or die "Cannot open $assfile: $!";
my %ashash = ();
my $ignore_good_no = 0;
my $ignore_poor_no = 0;
my $update_no = 0;
my $nonupdate_no = 0;
my $count = 0;
my %check = ();

<SRC>;

foreach my $line (<SRC>) {
	chomp $line;
	my @lines = split(/\t/, $line);
	
	if($lines[1] == 1 and $lines[2] == 0){
		if(not exists $check{$lines[0]}){
			$check{$lines[0]} = 0;
			$ignore_good_no ++;
		}
		next;
	}
	
	if($lines[2] == $lines[3]){
		if(not exists $check{$lines[0]}){
			$check{$lines[0]} = 0;
			$ignore_poor_no ++;
		}
		next;
	}
	
	my @go = split(/,/, $lines[5]);
	if($go[0] eq "NA"){
		$ashash{$lines[0]}{$lines[4]} = [];
	}else{
		$ashash{$lines[0]}{$lines[4]} = [@go];
	}
}

foreach my $group (sort(keys %ashash)){
	my @go_list = ();
	my @longest = ();
	my @intersect = ();
	my $go_comp = 0;
	my $mark = 0;
	my @ATs = grep(/^AT/, sort(keys %{$ashash{$group}}));
	
	my @genes = (@ATs)? @ATs : sort(keys %{$ashash{$group}});
	
	foreach my $gene (@genes){
		if(@{$ashash{$group}{$gene}}){
			push @go_list, \@{$ashash{$group}{$gene}};
			@longest = @{$ashash{$group}{$gene}} if (scalar(@longest) < scalar(@{$ashash{$group}{$gene}}));
		}
	}
	
	if(scalar(@go_list) > 1){
		$go_comp = List::Compare->new(@go_list);
		@intersect = $go_comp->get_intersection;
	}
				
	foreach my $short (@go_list) {
		$go_comp = List::Compare->new($short, \@longest);
		if($go_comp->is_LsubsetR){
			$mark = 1;
		}else{
			$mark = 2;
		}
	}
	#print "$mark\n";
	
	if($mark){
		$update_no ++;
		foreach my $gene (sort(keys %{$ashash{$group}})) {
			my ($prefix) = $gene =~ /^([A-Za-z]+)/;			## get prefix
			my $sp = &prefix2species($prefix);				## get sp name
			if($mark == 1){
				$bighash{$sp}{$gene}{"GO"} = [@longest];
			}elsif($mark == 2){
				$bighash{$sp}{$gene}{"GO"} = [@intersect];
			}
		}
	}else{
		$nonupdate_no ++;
	}
}

print "Ortholog groups ignored due to good annotations: ", $ignore_good_no, "\n";
print "Ortholog groups ignored due to most genes unannotated: ", $ignore_poor_no, "\n";
print "Ortholog groups ignored due to unsolved inconsistency: ", $nonupdate_no, "\n";
print "Ortholog groups updated: ", $update_no, "\n";

close(SRC);

###########################################################################################
###  Generate edited annotation file with updated annot info from other species which	###
###  have better annot.																	###
###########################################################################################

print "******************** Generate edited annotation file *******************\n";

@sp = ("ath", "ptr", "gma", "pvu");
for (my $i=0; $i<scalar(@sp); $i++) {											## for each specise
	open(ANT, ">".$tgtannot[$i]) or die "Cannot write ".$tgtannot[$i]." :$!";   ## open new annot file to write
	
	foreach my $key (sort(keys $bighash{$sp[$i]})) {							## for each gene
		my %gannot = %{$bighash{$sp[$i]}{$key}};					
		my $ID = (@{$gannot{"ID"}}) ? join(",",@{$gannot{"ID"}}) : "";					## Phytozome internal transcript ID
		my $TName = (@{$gannot{"TName"}}) ? join(",",@{$gannot{"TName"}}) : "";			## Phytozome transcript name
		my $PName = (@{$gannot{"PName"}}) ? join(",",@{$gannot{"PName"}}) : "";			## Phytozome protein name
		my $PFAM = (@{$gannot{"PFAM"}}) ? join(",",@{$gannot{"PFAM"}}) : "";			## PFAM
		my $Panther = (@{$gannot{"Panther"}}) ? join(",",@{$gannot{"Panther"}}) : "";	## Panther
		my $KOG = (@{$gannot{"KOG"}}) ? join(",",@{$gannot{"KOG"}}) : "";				## KOG
		my $EC = (@{$gannot{"EC"}}) ? join(",",@{$gannot{"EC"}}) : "";					## KEGG ec
		my $KEGG = (@{$gannot{"KEGG"}}) ? join(",",@{$gannot{"KEGG"}}) : "";			## KEGG Orthology
		my $GO = (@{$gannot{"GO"}}) ? join(",",@{$gannot{"GO"}}) : "";					## Gene Ontology terms
		my $HitName = (@{$gannot{"HitName"}}) ? join(",",@{$gannot{"HitName"}}) : "";	## best arabidopsis TAIR10 hit name
		my $Symbol = (@{$gannot{"Symbol"}}) ? join(",",@{$gannot{"Symbol"}}) : "";		## best arabidopsis TAIR10 hit symbol
		my $HD = $gannot{"HitDefline"} ? $gannot{"HitDefline"} : "";	## best arabidopsis TAIR10 hit defline
		print ANT "$ID\t$key\t$TName\t$PName\t$PFAM\t$Panther\t$KOG\t$EC\t$KEGG\t$GO\t$HitName\t$Symbol\t$HD\n";
	} ## foreach

	close(ANT);
} ## for

print "******************** Done *******************\n\n";

system("ln -sf ./Athaliana_167_TAIR10.annotation_edit_V3.0.txt ./arabidopsis/00.data/annotation_info_latest.tab");
system("ln -sf ./Ptrichocarpa_210_v3.0.annotation_edit_V3.0.txt ./populus/00.data/annotation_info_latest.tab");
system("ln -sf ./Gmax_189_annotation_edit_V3.0.txt ./soybean/00.data/annotation_info_latest.tab");
system("ln -sf ./Pvulgaris_218_annotation_edit_V3.0.txt ./common.bean/00.data/annotation_info_latest.tab");

###########################################################################################
###  Defined functions.																	###
###########################################################################################

sub prefix2species{													## gene prefix to sp name
	my $pf = shift;
	$pf = lc($pf);
	return "ath" if $pf =~ /^at/;
	return "gma" if $pf eq "glyma";
	return "ptr" if $pf eq "potri";
	return "pvu" if $pf eq "phvul";
}