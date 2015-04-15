#!/usr/bin/perl -w

### GLOBAL VARS ###
###################
my $infile = $ARGV[0];
### MAIN CODE ###
#################
open(FIN,"$infile") || die "Can't open $infile: $!\n";

while(<FIN>){
    
    chomp; my $line = $_; 
    my ($chr, $cds_start, $cds_stop, $ID, $gene, $exon_num, $exon_start, $exon_stop) = split(/\s/, $line);
    unless($gene){$gene = " " ;}
    $chr=~s/^chr//;
    my ($id, $version) = split(/\./, $ID); #for ccds data
    my @chrarray = split(/\_/, $chr);
    my $chrNum = $chrarray[0];
    my $chrcount = scalar(@chrarray);
    if ($chrcount == 1){  #if need the halytype like chr6_hap, need to screen out
	my @exon_startarray = split(/\,/, $exon_start);
	my @exon_stoparray = split(/\,/, $exon_stop);
	my $outfile1 = "${infile}-exon.bed";
	my $outfile2 = "${infile}-Coding_region.bed";
	open(OUT1, ">>$outfile1") || die "Can't open $outfile1: $!\n";
	#open(OUT2, ">>$outfile2") || die "Can't open $outfile2: $!\n";
	#print OUT2 "$chrNum\t$cds_start\t$cds_stop\t${id}\n" if($cds_start != $cds_stop);
	#print OUT2 "${chrNum}\t$cds_start\t$cds_stop\t${gene}\n";
	for (my $j=0; $j<$exon_num; $j++) { #for gene
	#for (my $j=0; $j<$exon_num; $j++){ #for ;
	    next if ($exon_startarray[$j] > $cds_stop);
	    next if ($exon_stoparray[$j]< $cds_start);
	    
	    
	    if($cds_start != $cds_stop) {
		if ($exon_startarray[$j] >= $cds_start && $exon_stoparray[$j] <= $cds_stop) {
		    print OUT1 "$chrNum\t$exon_startarray[$j]\t$exon_stoparray[$j]\t$gene\t${id}_exon_${j}\n";
		}elsif($exon_startarray[$j] < $cds_start && $exon_stoparray[$j] <= $cds_stop) {
		    print OUT1 "$chrNum\t$cds_start\t$exon_stoparray[$j]\t$gene\t${id}_exon_${j}\n";
		}elsif($exon_startarray[$j] >= $cds_start && $exon_stoparray[$j] > $cds_stop) {
		    print OUT1 "$chrNum\t$exon_startarray[$j]\t$cds_stop\t$gene\t${id}_exon_${j}\n";
		}elsif($exon_startarray[$j] < $cds_start && $exon_stoparray[$j] > $cds_stop) {
                    print OUT1 "$chrNum\t$cds_start\t$cds_stop\t$gene\t${id}_exon_${j}\n";
                }
	    }else{
		print OUT1 "$chrNum\t$exon_startarray[$j]\t$exon_stoparray[$j]\t$gene\t${id}_exon_${j}\n";
	    }
	    
	}
    }
    
}


close(FIN);
close(OUT1);
#close(OUT2);
