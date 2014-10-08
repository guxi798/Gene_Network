#!/usr/bin/perl -w
use strict;

my $sp = shift @ARGV;
my $annotDB = shift @ARGV;
my $gmtFile = "01.GeneSet.".$annotDB;

my $srcfolder = "./".$sp."/03.result/08.GSEA";
opendir(SRC, $srcfolder);
my @files = readdir(SRC);
my @txt = sort(grep(/02.ExpressionSet.datExpr/, @files));
my @rnk = sort(grep(/04.RankedGeneList.Eigengene/, @files));
my @gmt = sort(grep(/$gmtFile/, @files));

my $queue = "rcc-mc-30d";
my @para = ();
push @para,("-cp", "01.script/gsea2.jar", "-Xmx5000m", "xtools.gsea.GseaPreranked");

push @para,("-nperm", "1000");
push @para,("-collapse", "false");
#push @para,("-permute", "gene_set");
push @para,("-scoring_scheme", "classic");
#push @para,("-metric", "Pearson");
#push @para,("-sort", "real");
#push @para,("-order", "descending");
push @para,("-set_max", "500");
push @para,("-set_min", "5");
push @para,("-norm", "meandiv");
push @para,("-make_sets", "true");
#push @para,("-median", "false");
#push @para,("-num", "100");
push @para,("-plot_top_x", "20");
#push @para,("-rnd_seed", "timestamp");
push @para,("-rnd_seed", "149");
#push @para,("-save_rnd_lists", "false");
push @para,("-zip_report", "false");
push @para,("-out", $srcfolder."/preranked");
push @para,("-gui", "false");

system("mkdir -p $srcfolder/preranked");
system("rm -r $srcfolder/preranked/*");

for(my $i=0; $i<scalar(@txt); $i++){
    my ($module) = $txt[$i] =~ /02\.ExpressionSet\.datExpr\.([a-zA-Z]+)\.txt/;
    
    my $tgtfile = "GSEA.folderParameters.preranked.$sp.$module.sh";
    open(TGT, ">$tgtfile");
    print TGT "#!/bin/bash\n\n";
    print TGT "time java ", join(" ", @para), " -rpt_label ", $module, " -rnk $srcfolder/$rnk[$i] -gmx $srcfolder/$gmt[0]\n";
    close(TGT);

    system("chmod u+x $tgtfile");
    system("qsub -q $queue -e 01.script/e.GSEA.preranked.$sp.$module -o 01.script/o.GSEA.preranked.$sp.$module ./$tgtfile");
    system("rm $tgtfile");
}

close(SRC);





