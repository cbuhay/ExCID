#!/usr/bin/perl -w
use strict;


my $annotated_index = $ARGV[0];
my $HGNC = $ARGV[1];


open(my $fh,"<$annotated_index") or die $!;


while (my $line = <$fh>) {
    
    chomp($line);
    
    my ($chr,$start,$Stop,$transcript_ID) = split("\t",$line);
    
    my @grep = `grep -w "$transcript_ID" $HGNC `;
    
    if (scalar(@grep) == 1) {
        my @tmp = split("\t",$grep[0]);
        my $gene_name = $tmp[0];
        print "$chr\t$start\t$Stop\t$transcript_ID\t$gene_name\n";
    }else{
        #my @tmp = split("-",$transcript_ID);
        #my $gene_name = "MIR";
        #for(my $i = 2; $i < scalar(@tmp); $i++){
        #    $gene_name.=uc($tmp[$i]);
        #}
        print "$chr\t$start\t$Stop\t$transcript_ID\t.\n";
        #print STDERR "$chr\t$start\t$Stop\t$gene_name\t$transcript_ID\n";
    }
    
}

close($fh);