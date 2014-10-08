###########################################################
###	@Author		Xi Gu									###
###	@Time		Aug 05, 2014							###
###	@Function	Assess KEGG Annotations					###
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

my @sourcefiles = qw(./arabidopsis/00.data/Athaliana_167_TAIR10.annotation_info.txt
				   ./populus/00.data/Ptrichocarpa_210_v3.0.annotation_info.txt
				   ./soybean/00.data/Gmax_189_annotation_info.txt
				   ./common.bean/00.data/Pvulgaris_218_annotation_info.txt);
my @tgtannot = qw(./arabidopsis/00.data/Athaliana_167_TAIR10.annotation_edit_V1.0.txt
				   ./populus/00.data/Ptrichocarpa_210_v3.0.annotation_edit_V1.0.txt
				   ./soybean/00.data/Gmax_189_annotation_edit_V1.0.txt
				   ./common.bean/00.data/Pvulgaris_218_annotation_edit_V1.0.txt);
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
		} # if
	} # foreach
	close(SRC);

	my $sp = shift(@sp);								## species name
	$bighash{$sp} = {%sphash};							## key is sp name
} # foreach

print "******************** Done *******************\n\n";

###########################################################################################
###  Provide location of ortholog group file and output files to write assess result.	###
###  Here we are only interested in KEGG and GO databases.								###
###########################################################################################

##### open file handles for ortholog group file, KEGG assess result file and GO assess result file. #####
my $ogfile = "./common.bean/00.data/ortholog_groups.tab";
my $tgtfile1 = "./common.bean/03.result/09.annotAssessment/01.annotAssessment.KEGG.tab";

open(SRC, "$ogfile");
open(KEG, ">$tgtfile1");

##### print title line to the two assessment files.	#####
print KEG "mcl_group\t", "JC\t","unannotated_No\ttotal_No\t","unnatated_ath\ttotal_ath\t",
		"unnatated_ptr\ttotal_ptr\t", "unnatated_gma\ttotal_gma\t", "unnatated_pvu\ttotal_pvu\n";

###########################################################################################
###  Start assessing the annotation. The basic criterion is to evaluate the extent to	###
###  which annotations of different species ortholog genes agree with each other. 		###
###  For Kegg, as annotations are generated from Ath, as long as the orthologs have		###
###  Kegg assigned, they must agree with each other. The only error will be some ortho-	###
###  -logs missed Kegg annotation. In this case, we will assign those genes with KEGG	###
###  KO terms from other annotated orthologs in the same group.							###
###  For GO, different genome can be annotated by different people, thus there can have ###
###  discrepancy among ortholog genes. In this case, as it cannot simply decide which   ###
###  annotation is correct, we will first pick out and print the inconsistent OGs for   ###
###  manual check.																		###
###########################################################################################

print "******************** Assess Annotation *******************\n";

my $i = 0;
while( my $line = <SRC>) {								## each line in OG file
	#print $i++, "\n";
	chomp $line;
	my @lines = split(/\t/, $line);
	my $group = shift @lines;							## group name

	my @kegg_null = ();									## count No of genes without KEGG annot
	my @kegg_store = ();								## store array to get intersection and union

	my @kegg_annot = ();

	my @kegg_inter = ();								## annot interaction among genes within the same OG
	my @kegg_union = ();								## annot union among genes within the same OG

	my %kegg_sp_null = ();								## count No of genes without KEGG annot in each sp
	my %full = ();										## for counting total no of genes and genes in each sp

##### initialize hash for each sp as the key #####
	foreach my $sp (@sp) {								
		$full{$sp} = ();
		$kegg_sp_null{$sp} = ();
	} ## foreach

##### initial object to obtain intersection/union and store genes with no annotations ######
	foreach my $gene (@lines) {							## for each ortholog
		my ($prefix) = $gene =~ /^([A-Za-z]+)/;			## get prefix
		my $sp = &prefix2species($prefix);				## get sp name
		push @{$full{$sp}}, $gene;

		if(@{$bighash{$sp}{$gene}{"KEGG"}}){						## if there is KEGG annot available for this gene
			push @kegg_store, \@{$bighash{$sp}{$gene}{"KEGG"}};		## store array addresses to construct comp object
			if(scalar(@kegg_annot) < scalar(@{$bighash{$sp}{$gene}{"KEGG"}})){
				@kegg_annot = @{$bighash{$sp}{$gene}{"KEGG"}};
			} ## if
		}else{														## no KEGG annot available for this gene
			push @kegg_null, $gene;									## put into null array
			push @{$kegg_sp_null{$sp}}, $gene;						## put into sp hash of null array
		} ## if
	} ## foreach

##### obtain intersection and union ######
	if(@kegg_annot){
		foreach my $gene (@lines) {
			my ($prefix) = $gene =~ /^([A-Za-z]+)/;			## get prefix
			my $sp = &prefix2species($prefix);				## get sp name
			my $lr = List::Compare->new(\@kegg_annot, \@{$bighash{$sp}{$gene}{"KEGG"}}); 

			if(!$lr->is_LeqvlntR){
				@{$bighash{$sp}{$gene}{"KEGG"}} = @kegg_annot;
			} ## if
		} ## foreach
	} ## if

##### compute intersection and union #####
	#print join(",", @kegg_store), "\n";
	my $kegg_J = 0;
	if (scalar(@kegg_store) > 1) {
		my $kegg_comp = List::Compare->new(@kegg_store);		## construct new object
		@kegg_inter = $kegg_comp->get_intersection;						## get interesction
		@kegg_union = $kegg_comp->get_union;							## get union
		$kegg_J = scalar(@kegg_inter)/scalar(@kegg_union);
	}

##### print results into KEGG output file #####
	print KEG "$group\t", $kegg_J, "\t", scalar(@kegg_null), "\t", scalar(@lines)+1, "\t";
	@sp = ("ath", "ptr", "gma", "pvu");
	foreach my $sp (@sp) {
		if($kegg_sp_null{$sp}){
			print KEG scalar @{$kegg_sp_null{$sp}}, "\t", scalar @{$full{$sp}}, "\t";
		}else{
			print KEG 0, "\t";
			if($full{$sp}){
				print KEG scalar @{$full{$sp}}, "\t";
			}else{
				print KEG 0, "\t";
			} ## if
		} ## if
	} ## foreach
	print KEG "\n";
} # while

##### close file handles #####
close(SRC);
close(KEG);

print "******************** Done *******************\n\n";

###########################################################################################
###  Generate edited annotation file with updated annot info from other species which	###
###  have better annot.																	###
###########################################################################################

print "******************** Generate edited annotation file *******************\n";

for (my $i=0; $i<scalar(@sp); $i++) {											## for each specise
	open(ANT, ">".$tgtannot[$i]) or die "Cannot write ".$tgtannot[$i]." :$!";   ## open new annot file to write
	
	foreach my $key (sort(keys $bighash{$sp[$i]})){							## for each gene
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

system("ln -sf ./Athaliana_167_TAIR10.annotation_edit_V1.0.txt ./arabidopsis/00.data/annotation_info_latest.tab");
system("ln -sf ./Ptrichocarpa_210_v3.0.annotation_edit_V1.0.txt ./populus/00.data/annotation_info_latest.tab");
system("ln -sf ./Gmax_189_annotation_edit_V1.0.txt ./soybean/00.data/annotation_info_latest.tab");
system("ln -sf ./Pvulgaris_218_annotation_edit_V1.0.txt ./common.bean/00.data/annotation_info_latest.tab");

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