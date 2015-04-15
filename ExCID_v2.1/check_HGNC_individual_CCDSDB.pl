#!/usr/bin/perl -w
use strict;


my $annotated_index = $ARGV[0];
my $HGNC = $ARGV[1];


open(my $fh,"<$annotated_index") or die $!;


while (my $line = <$fh>) {
    
    chomp($line);
    
    my ($chr,$start,$Stop,$ID) = split("\t",$line);
    my @transcript_ID_split = split("_exon_",$ID);
    my $transcript_ID = $transcript_ID_split[0];
    
    my @grep = `grep -w "$transcript_ID" $HGNC `;
    
    if (scalar(@grep) == 1) {
        my @tmp = split("\t",$grep[0]);
        my $gene_name = $tmp[0];
        print "$chr\t$start\t$Stop\t$ID\t$gene_name\n";
    }else{
        print "$chr\t$start\t$Stop\t$ID\t.\n";
        #print STDERR "$chr\t$start\t$Stop\t$ID\t.\n";
    }
    
}

close($fh);
