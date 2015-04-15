#!/usr/bin/perl -w
# This script takes in the annotated BED file (eg the VCrome+PKv2 fully annotated) and checks
# the gene names and other gene names columns with the HGNC database to add the Approved names in 
# the gene names column and put other names in the other gene names column.
#

use strict;


my $annotated_index = $ARGV[0];
my $HGNC = $ARGV[1];


open(my $fh,"<$annotated_index") or die $!;


while (my $line = <$fh>) {
    
    chomp($line);
    
    $line=~s/\t-1\t/\t.\t/;
    my $not_anno = 0;
    my @line_tmp = split("\t",$line);
    my $target = $line_tmp[0]."\t".$line_tmp[1]."\t".$line_tmp[2];
    my $gene_name="";
    my $prev_name="";
    my $synonyms_names = "";
    my $refseq_IDs = $line_tmp[11];
    my $CCDS_IDs = $line_tmp[12];
    my $VEGA_IDs = $line_tmp[13];
    my $miRNA_IDs = $line_tmp[14];
    my $rest = join("\t", @line_tmp[15..(scalar(@line_tmp)-1)]);
    
    
    ## Gene_name
    
    if ($line_tmp[3] ne "." ) {
        $gene_name =$line_tmp[3].";";
        my $check = $line_tmp[4].";";
        if(index($synonyms_names,$check) == -1 && index($gene_name,$check) == -1) {
            $synonyms_names = $line_tmp[4].";";   ## First time Syn Name is given a value.
        }
        $check = $line_tmp[5].";";
        if(index($synonyms_names,$check) == -1 && index($gene_name,$check) == -1) {
            $synonyms_names .= $line_tmp[5].";";
        }
        $check = $line_tmp[6].";";
        if(index($synonyms_names,$check) == -1 && index($gene_name,$check) == -1) {
            $synonyms_names .= $line_tmp[6].";";
        }
    }elsif($line_tmp[4] ne "." ){
        $gene_name =$line_tmp[4].";";
        my $check = $line_tmp[5].";";
        if(index($synonyms_names,$check) == -1 && index($gene_name,$check) == -1) {
            $synonyms_names = $line_tmp[5].";"; ## First time Syn Name is given a value.
        }
        $check = $line_tmp[6].";";
        if(index($synonyms_names,$check) == -1 && index($gene_name,$check) == -1) {
            $synonyms_names .= $line_tmp[6].";";
        }
    }elsif($line_tmp[5] ne "." ){
        $gene_name =$line_tmp[5].";";
        my $check = $line_tmp[6].";";
        if(index($synonyms_names,$check) == -1 && index($gene_name,$check) == -1) {
            $synonyms_names = $line_tmp[6].";"; ## First time Syn Name is given a value.
        }
    }elsif($line_tmp[6] ne "." ){
        $gene_name =$line_tmp[6].";";
    }
    
    ##
    
    ## Other_Names
    
    if ($line_tmp[7] ne "."){
        my @line_tmp_split = split(";",$line_tmp[7]);
        foreach my $tmp (@line_tmp_split){
            my $check = $tmp.";";
            if (index($gene_name,$check) == -1 && index($synonyms_names,$check) == -1 && index($prev_name,$check) == -1) {
                $prev_name = $line_tmp[7].";";
            }
            
        } 
    }
    
    if ($line_tmp[8] ne "."){
        my @line_tmp_split = split(";",$line_tmp[8]);
        foreach my $tmp (@line_tmp_split){
            my $check = $tmp.";";
            if (index($gene_name,$check) == -1 && index($synonyms_names,$check) == -1 && index($prev_name,$check) == -1) {
                $prev_name = $line_tmp[8].";";
            }
            
        } 
    }
    
    if ($line_tmp[9] ne "."){
        my @line_tmp_split = split(";",$line_tmp[9]);
        foreach my $tmp (@line_tmp_split){
            my $check = $tmp.";";
            if (index($gene_name,$check) == -1 && index($synonyms_names,$check) == -1 && index($prev_name,$check) == -1) {
                $prev_name = $line_tmp[9].";";
            }
            
        } 
    }
    
    if ($line_tmp[10] ne "."){
        my @line_tmp_split = split(";",$line_tmp[10]);
        foreach my $tmp (@line_tmp_split){
            my $check = $tmp.";";
            if (index($gene_name,$check) == -1 && index($synonyms_names,$check) == -1 && index($prev_name,$check) == -1) {
                $prev_name = $line_tmp[10].";";
            }
            
        } 
    }
    
    $gene_name=~ s/;$//;
    $prev_name=~ s/;$//;
    $synonyms_names=~ s/;$//;
    $gene_name=~ s/^\.;//;
    $prev_name=~ s/^\.;//;
    $synonyms_names=~ s/^\.;//;
    $gene_name=~ s/;\.$//;
    $prev_name=~ s/;\.$//;
    $synonyms_names=~ s/;\.$//;
    
    chomp($gene_name);
    chomp($prev_name);
    chomp($synonyms_names);
    
    if (length($gene_name) == 0) {
        $gene_name = ".";
    }
    if (length($prev_name) == 0) {
        $prev_name = ".";
    }
    if (length($synonyms_names) == 0) {
        $synonyms_names = ".";
    }
    
    print "$target\t$gene_name\t$prev_name\t$synonyms_names\t$refseq_IDs\t$CCDS_IDs\t$VEGA_IDs\t$miRNA_IDs\t$rest\n";
 
}

close($fh);
