###########################################################
###	@Author		Xi Gu									###
###	@Time		July 21, 2014							###
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
my @uniprot = qw(./arabidopsis/00.data/Uniprot_Arabidopsis.tab
				 ./populus/00.data/Uniprot_Populus.tab
				 ./soybean/00.data/Uniprot_Glycine.tab);
my @sp = ("ath", "ptr", "gma", "pvu");					## species name
my $ptrconver = "./populus/00.data/Ptrichocarpa_210_v3.0.synonym.txt";

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

print "******************** Add UniProt annotations to the existing structure *******************\n";
@sp = ("ath", "ptr", "gma", "pvu");
foreach my $file (@uniprot) {						## for each species
	my $sp = shift(@sp);								## species name
	my %ptrconver = ();
	open(SRC, "$file") or die "Cannot open $file: $!";

	<SRC>;
	my $mark = 0;
	$mark = 1 if($sp eq "gma");
	$mark = 2 if($sp eq "ptr");
	
	if($mark == 2){
		open(PTR, "$ptrconver") or die "Cannot open $ptrconver: $!";
		
		foreach my $line (<PTR>) {
			chomp $line;
			my @lines = split(/\t/, $line);
			my @newids = split(/\./, $lines[0]);
			pop @newids;
			my $newid = join(".", @newids);
			
			$ptrconver{$lines[1]} = $newid if(not exists $ptrconver{$lines[1]});
		}
		
		close(PTR)
	}
	
	foreach my $line (<SRC>) {							## for each gene entry
		chomp $line;
		my @lines = split(/\t/, $line);
		my @ids = split(/\./, $lines[7]);
		pop @ids;
		my $id = join(".", @ids);
		
		if($mark == 1){
			$id = lc($id);
			$id =~ s/glyma/Glyma/;
		}
		if($mark == 2){
			$id = $ptrconver{$id} if(exists $ptrconver{$id});
		}
		if (exists $bighash{$sp}{$id}) {			## key is gene name
			#print "$sp\t$id\t", $lines[3], "\n";
			$bighash{$sp}{$id}{"UniProt"} = $lines[3];
		}
	}
	close(SRC);
}

###########################################################################################
###  Provide location of ortholog group file and output files to write assess result.	###
###  Here we are only interested in KEGG and GO databases.								###
###########################################################################################

##### open file handles for ortholog group file, KEGG assess result file and GO assess result file. #####
my $ogfile = "./common.bean/00.data/ortholog_groups.tab";
my $tgtfile2 = "./common.bean/03.result/09.annotAssessment/02.annotAssessment.GO.tab";
my $tgtfile3 = "./common.bean/03.result/09.annotAssessment/03.annotAssessment.GO.detail.tab";

open(SRC, "$ogfile");
open(GOT, ">$tgtfile2");
open(GOD, ">$tgtfile3");

##### print title line to the two assessment files.	#####
print GOT "mcl_group\t","JC\t","unannotated_No\ttotal_No\t","unnatated_ath\ttotal_ath\t",
		"unnatated_ptr\ttotal_ptr\t", "unnatated_gma\ttotal_gma\t", "unnatated_pvu\ttotal_pvu\n";
print GOD "mcl_group\t","JC\t","unannotated_No\ttotal_No\t", "GO_terms\t", "Phytozome_description\tUniprot_description\n";

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

	my @go_null = ();									## count No of genes without GO annot
	my @go_store = ();

	my @go_inter = ();
	my @go_union = ();

	my %go_sp_null = ();								## count No of genes without GO annot in each sp
	my %full = ();										## for counting total no of genes and genes in each sp

##### initialize hash for each sp as the key #####
	foreach my $sp (@sp) {								
		$full{$sp} = ();
		$go_sp_null{$sp} = ();
	} ## foreach

##### initial object to obtain intersection/union and store genes with no annotations ######
	foreach my $gene (@lines) {							## for each ortholog
		my ($prefix) = $gene =~ /^([A-Za-z]+)/;			## get prefix
		my $sp = &prefix2species($prefix);				## get sp name
		push @{$full{$sp}}, $gene;

		if(@{$bighash{$sp}{$gene}{"GO"}}){							## if there is GO annot available for this gene
			push @go_store, \@{$bighash{$sp}{$gene}{"GO"}};			## store array addresses to construct comp object
		}else{														## no GO annot available for this gene
			push @go_null, $gene;									## put into null array
			push @{$go_sp_null{$sp}}, $gene;						## put into sp hash of null array
		}  ## if
	} ## foreach

##### compute intersection and union #####
	my $go_J = 0;
	if (scalar(@go_store) > 1) {
		my $go_comp = List::Compare->new({lists => [@go_store]});	
		@go_inter = $go_comp->get_intersection;
		@go_union = $go_comp->get_union;
		$go_J = scalar(@go_inter)/scalar(@go_union);
	}

##### print results into GO output file #####
	print GOT "$group\t", $go_J,"\t", scalar(@go_null), "\t", scalar(@lines), "\t";
	@sp = ("ath", "ptr", "gma", "pvu");
	foreach my $sp (@sp) {
		if($go_sp_null{$sp}){
			print GOT scalar @{$go_sp_null{$sp}}, "\t", scalar @{$full{$sp}}, "\t";
		}else{
			print GOT 0, "\t";
			if($full{$sp}){
				print GOT scalar @{$full{$sp}}, "\t";
			}else{
				print GOT 0, "\t";
			} ## if
		} ## if
	} ## foreach
	print GOT "\n";

##### print results into GO output file which contains specific GO assignment for each gene in each OG #####
	foreach my $gene (@lines) {
		print GOD "$group\t", $go_J,"\t", scalar(@go_null), "\t", scalar(@lines), "\t";
		my ($prefix) = $gene =~ /^([A-Za-z]+)/;			## get prefix
		my $sp = &prefix2species($prefix);				## get sp name
		print GOD "$gene\t";
		(@{$bighash{$sp}{$gene}{"GO"}}) ? print GOD join(",", @{$bighash{$sp}{$gene}{"GO"}}) : print GOD "NA";
		print GOD "\t";
		($bighash{$sp}{$gene}{"HitDefline"}) ? print GOD $bighash{$sp}{$gene}{"HitDefline"} : print GOD "NA";
		if($sp eq "pvu"){print GOD "\n"; next;}
		print GOD "\t";
		($bighash{$sp}{$gene}{"UniProt"}) ? print GOD $bighash{$sp}{$gene}{"UniProt"} : print GOD "NA";
		print GOD "\n";
	}
} ## while

##### close file handles #####
close(SRC);
close(GOT);
close(GOD);

print "******************** Done *******************\n\n";

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