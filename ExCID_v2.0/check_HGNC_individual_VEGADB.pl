#!/usr/bin/perl -w
use strict;


my $annotated_index = $ARGV[0];
my $HGNC = $ARGV[1];
my $VEGA_HGNC_names = $ARGV[2];

open(my $fh,"<$annotated_index") or die $!;


while (my $line = <$fh>) {
    
    chomp($line);
    
    my ($chr,$start,$Stop,$gene,$ID) = split("\t",$line);
    my @transcript_ID_tmp_split = split("_exon_",$ID);
    my $transcript_ID = $transcript_ID_tmp_split[0];
    
    my @grep2 = `grep -w "$transcript_ID" $VEGA_HGNC_names `;
    if (scalar(@grep2) > 0) {
        my @grep_tmp = split("\t",$grep2[0]);
        chomp($grep_tmp[0]);
        $transcript_ID= $grep_tmp[0];  
    }
    
    my @grep = `grep -w "$transcript_ID" $HGNC `;
    
    if (scalar(@grep) == 1) {
        my @tmp = split("\t",$grep[0]);
        my $gene_name = $tmp[0];
        print "$chr\t$start\t$Stop\t$gene\t$ID\t$gene_name\n";
    }else{
        print "$chr\t$start\t$Stop\t$gene\t$ID\t.\n";
       # print STDERR "$chr\t$start\t$Stop\t$gene\t$ID\t.\n";
    }
    
}

close($fh);
