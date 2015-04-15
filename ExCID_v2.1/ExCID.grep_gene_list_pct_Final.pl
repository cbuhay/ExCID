#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Path;
use Time::localtime;
use Fcntl qw(:flock);
use File::Basename;


unless (scalar @ARGV > 0){USAGE(); exit; }

my %opt;
my $script_dir;
my $script_dir_tmp = $0;
$script_dir_tmp =~m/^.+\//;
$script_dir=$&;

my @tmp = `grep "DataBaseDir=" $script_dir/Config.txt`;
$tmp[0]=~s/DataBaseDir=//;
chomp($tmp[0]);
my $database_dir = $tmp[0];

my $database1 = "$database_dir/refGene.txt.bed-exon_HGNC.bed";
my $database2 = "$database_dir/ccdsGene.txt.bed-exon_HGNC.bed";
my $database3 = "$database_dir/vegaGene.txt.bed-exon_HGNC.bed";
my $database4 = "$database_dir/miRBASE_r20_HGNC.bed";
my $combined='';
my $HGMD_db = "$database_dir/HGMD_2014_v4.bed";


### OPTIONS ###
my $check_HGMD;
GetOptions('i:s' => \$opt{i}, 'l:s' => \$opt{l}, 'list:s' => \$opt{g}, 'db:s' => \$opt{db}, 'checkHGMD' => \$check_HGMD) || &USAGE;
my $bed = $opt{i} || "null";
my $low_cov_file = $opt{l} || "null";
my $gene_list;
$gene_list = $opt{g} if($opt{g});
my $targeted_database;
$targeted_database = $opt{db} if($opt{db}) ;
my $pct_file;
my @genes_lest;


if (!$gene_list && !$targeted_database) {
    print STDERR "Please provide either list of genes (HGNC Symbols) or the targeted Gene database. Check the documentation for database format\n";
    exit;
}

if ($bed eq "null") {
    print STDERR "Please provide Targeted Bedfile.\n";
    exit;
}

if($low_cov_file eq "null"){
    print STDERR "Please provide inadequately covered bases in bed file format.\n";
    exit;
}



if($targeted_database && -e $targeted_database && !$gene_list){
    print STDERR  "Targeted Data base provided. Obtaining Gene percent coverages.\n";
    $pct_file=genes_check_db($low_cov_file, $bed, $targeted_database);
    if($check_HGMD) {
        HGMDcheck($pct_file,@genes_lest) ;
    }else{
        averageGene($pct_file,@genes_lest) ;
    }
    exit 0 ;
}elsif($targeted_database && -e $targeted_database && $gene_list && -e $gene_list){
    print STDERR  "Targeted Data base provided. Obtaining Gene percent coverages.\n";
    $pct_file=genes_check_db($low_cov_file, $bed, $targeted_database);
    open(my $glfh, "< $gene_list") or die $!;
    @genes_lest= <$glfh>;
    close($glfh);
    if($check_HGMD) {
        HGMDcheck($pct_file,@genes_lest) ;
    }else{
        averageGene($pct_file,@genes_lest) ;
    }
    exit 0 ;
}elsif($gene_list && -e $gene_list && !$targeted_database){
    open(my $glfh, "< $gene_list") or die $!;
    @genes_lest= <$glfh>;
    close($glfh);
}else{
    print STDERR "Some issue with the run. Please check the commands.\n";
}

my $output_dir;
my $output_dir_tmp = $low_cov_file;
$output_dir_tmp =~m/^.+\//;
$output_dir=$&;

$targeted_database = $output_dir."/".basename($gene_list)."-database.bed";
my $targeted_miRBASE = $output_dir."/".basename($gene_list)."-miRBASE.bed";
open(my $gldb, "> $targeted_database") or die $!;
open(my $gldbm, "> $targeted_miRBASE") or die $!;
my $mirfound = 0;
foreach my $gene (@genes_lest){
    chomp($gene);
    next unless(length($gene) > 0);
    my $found = 0;
    my @regions = `grep -w -P \"\t$gene\$\" $database1 `;
    if (scalar(@regions) > 0) {
        foreach my $reg (@regions){
            chomp($reg);
            my @tmp_split = split("\t",$reg);
            print $gldb "$tmp_split[0]\t$tmp_split[1]\t$tmp_split[2]\t$tmp_split[5]|$tmp_split[4]\n" ;    
        }
        $found = 1;
    }
    
    @regions = `grep -w -P \"\t$gene\$\" $database2 `;
    if (scalar(@regions) > 0) {
        foreach my $reg (@regions){
            chomp($reg);
            my @tmp_split = split("\t",$reg);
            print $gldb "$tmp_split[0]\t$tmp_split[1]\t$tmp_split[2]\t$tmp_split[4]|$tmp_split[3]\n" ;    
        }
        $found = 1;
    }
    
    @regions = `grep -w -P \"\t$gene\$\" $database3 `;
    if (scalar(@regions) > 0) {
        foreach my $reg (@regions){
            chomp($reg);
            my @tmp_split = split("\t",$reg);
            print $gldb "$tmp_split[0]\t$tmp_split[1]\t$tmp_split[2]\t$tmp_split[5]|$tmp_split[4]\n" ;    
        }
        $found = 1;
    }
    
    @regions = `grep -w -P \"\t$gene\$\" $database4 `;
    if (scalar(@regions) > 0) {
        foreach my $reg (@regions){
            chomp($reg);
            my @tmp_split = split("\t",$reg);
            print $gldbm "$tmp_split[0]\t$tmp_split[1]\t$tmp_split[2]\t$tmp_split[4]|$tmp_split[3]\n" ;    
        }
        $found = 1;
        $mirfound = 1;
    }
    
    if ($found != 1) {
        print STDERR "$gene is not a HGNC symbol. Please check the Gene or update the databases.\n";
    }
    
}

close($gldb);
close($gldbm);

$pct_file=genes_check_db($low_cov_file, $bed, $targeted_database);

open(my $glfh, "< $gene_list") or die $!;
@genes_lest= <$glfh>;
close($glfh);

if($check_HGMD) {
    HGMDcheck($pct_file,@genes_lest) ;
}else{
    averageGene($pct_file,@genes_lest) ;
}
    
if ($mirfound == 1) {
    system("$script_dir/ExCID.grep_gene_list_pct_Final-miRNA.pl -i $bed -l $low_cov_file -db $targeted_miRBASE ");
}



########## SUBROUTINES ###########
sub USAGE {
    print "\nUSAGE: $0 -i <Targeted BED file> -l <low cov file> -list <Genelist> -db <Gene database>\n";
    print "  -i:  Annotated bed file\n";
    print "  -l:  Low cov bed file\n";
    print "  -list:  List of genes. one per line\n";
    print "  -checkHGMD: Output only HGMD Transcripts if present or Average among all the transcripts.\n";
    print "  -db:  Database of interested Genes. Please see documentation for database format. \n";
    exit;
}

sub genes_check_db {
    
    my ($low_cov_file, $Bed_file, $db_list)  = @_;
    
    my $output_file_name = basename($low_cov_file);
    $output_file_name=~ s{\.[^.]+$}{}; # removes extension
    my $output_dir;
    my $output_dir_tmp = $low_cov_file;
    $output_dir_tmp =~m/^.+\//;
    $output_dir=$&;
    my $db_name = basename($db_list);
    my $out_file = $output_file_name;
    $out_file=~ s{.*/}{}; # remove path
    $out_file.="-$db_name.bed";
    my $nottargeted = $Bed_file;
    $nottargeted=~ s{.*/}{}; # remove path
    $nottargeted.="-notTrgtdin-$db_name-$output_file_name.bed";
    my $tmp = $out_file;
    $tmp=~ s{.*/}{}; # remove path
    $tmp=~s/\.bed$//;
    $combined = $nottargeted;
    $combined=~s/\-$output_file_name//;
    $combined=~s/\.bed$//;
    $combined.="-$tmp.bed";
    
    system("$script_dir/bin/bedtools subtract -a $db_list -b $Bed_file > $nottargeted");
    system("$script_dir/bin/bedtools intersect -a $db_list -b $low_cov_file > $out_file");
    system("cat $nottargeted $out_file | $script_dir/bin/bedtools sort -i  > $combined");
    system("$script_dir/bin/bedtools merge -i $combined > $combined-tmp");
    system("$script_dir/bin/bedtools intersect -a $db_list -b $combined-tmp > $combined");
    my $db_list_transcript_size = get_Transcirpt_size_db($low_cov_file,$db_list);
    my $out_FILE = get_Transcirpt_size($low_cov_file,$combined, $db_list_transcript_size);
    
    my $rm_file = "$combined-tmp";
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = $out_file;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = $nottargeted;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = "$combined";
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = $db_list_transcript_size;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file=~s/\.txt$/\_perExon\.txt/;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    
    
    $db_name=~ s{\.[^.]+$}{}; # removes extension
    system("mv $out_FILE $output_dir/$output_file_name-$db_name-pct.txt");
    $out_FILE=~s/\_pct\.txt$/\_perExon\_pct\.txt/;
    system("mv $out_FILE $output_dir/$output_file_name-$db_name-perExon_pct.txt");
    return "$output_dir/$output_file_name-$db_name-pct.txt";
}
    
sub get_Transcirpt_size{
    my ($low_cov_file,$combined,$db)=@_;
    
    my $tmp = $combined;
    $tmp=~ s{.*/}{}; # remove path
    $tmp=~s/\.bed$//;
    
    
    my $outfile_final = "$tmp"."_transcriptSIZE.txt";
    my $outfile_final_exon = "$tmp"."_transcriptSIZE_perExon.txt";
    my %sizehash;
    open(my $in,"$combined") || die "Can't open $combined: $!\n";
    open(my $out, ">$outfile_final") || die "Can't open $outfile_final: $!\n";
    open(my $outex, ">$outfile_final_exon") || die "Can't open $outfile_final_exon: $!\n";
    
    while(<$in>)
         {
             chomp; my $line = $_;
             next unless (length($line) != 0);
             my ($chr, $start, $stop, $gene) = split(/\s/, $line);
             my $exonsize = $stop - $start + 1;
             my @tranarray = split(/\;/, $gene);
             my $arraysize = scalar(@tranarray) - 1;
             for (my $j=0; $j<=$arraysize; $j++) {
                print $outex "$tranarray[$j]\t$exonsize\n";
                my @tmp = split(/\_/, $tranarray[$j]);
                my $unit = join("_",@tmp[0..(scalar(@tmp)-3)]);
                my $cds = $tmp[(scalar(@tmp)-2)];
                my $exon = $tmp[(scalar(@tmp)-1)];
                push @{$sizehash{$unit}}, $exonsize;
             }
         }
    
    foreach my $key ( keys %sizehash )
    {
        my $total_size = eval join '+', @{$sizehash{$key}};
        print $out "$key\t$total_size\n";
    }
    
    close($in);
    close($out);
    close($outex);
    my $outFile = get_pct($outfile_final,$db);
    my $outFile_perExon = get_pct_perExon($outfile_final,$db);
    return $outFile;
}

sub get_Transcirpt_size_db{
    my ($low_cov_file, $db_list)=@_;
    
    my $tmp = $low_cov_file;
    $tmp=~ s{\.[^.]+$}{}; # removes extension
    
    my $tmp1 = basename($db_list);
    $tmp.="-$tmp1";
    $tmp=~s/\.bed$//;
    
    my $outfile_final = "$tmp"."_transcriptSIZE.txt";
    my $outfile_final_exon = "$tmp"."_transcriptSIZE_perExon.txt";
    my %sizehash;
    open(my $in,"< $db_list") || die "Can't open $db_list: $!\n";
    open(my $out, ">$outfile_final") || die "Can't open $outfile_final: $!\n";
    open(my $outex, ">$outfile_final_exon") || die "Can't open $outfile_final_exon: $!\n";
    
    while(<$in>)
         {
             chomp; my $line = $_;
             unless(length($line) != 0){next;}
             my ($chr, $start, $stop, $gene) = split(/\s/, $line);
             my $exonsize = $stop - $start + 1;
             my @tranarray = split(/\;/, $gene);
             my $arraysize = scalar(@tranarray) - 1;
             for (my $j=0; $j<=$arraysize; $j++) {
                print $outex "$tranarray[$j]\t$exonsize\t$chr\t$start\t$stop\n";
                my @tmp = split(/\_/, $tranarray[$j]);
                my $unit = join("_",@tmp[0..(scalar(@tmp)-3)]);
                my $cds = $tmp[(scalar(@tmp)-2)];
                my $exon = $tmp[(scalar(@tmp)-1)];
                push @{$sizehash{$unit}{'exonsize'}}, $exonsize;
                $sizehash{$unit}{'chr'} = $chr;
             }
         }
    
    foreach my $key ( keys %sizehash )
    {
        my $total_size = eval join '+', @{$sizehash{$key}{'exonsize'}};
        print $out "$key\t$total_size\t".$sizehash{$key}{'chr'}."\t".scalar(@{$sizehash{$key}{'exonsize'}})."\n";
        my @tmp= split(/\|/,$key);
        push @genes_lest,$tmp[0];
    }
    
    close($in);
    close($out);
    close($outex);
    return $outfile_final; 
}

sub get_pct{
    
    my ($infile,$control) =@_;
    
    unless (-e $control) {print STDERR "$control does not exist\n"; exit;}
    unless (-e $infile) {print STDERR "$infile does not exist\n"; exit;}
    
    my $outfile = $infile;
    $outfile=~s/\_transcriptSIZE\.txt$/\_pct\.txt/;
    my %control;
    my %file_data;
    open(my $fh,"<$control") or die $!;
    
    while (my $line = <$fh>) {
        chomp($line);
        my @data1 = split("\t",$line);
        #print STDERR $data1[0]."\n";
        my @data = split(/\|/,$data1[0]);
        #print STDERR $data[0]."\n";
        if (scalar @data1 ==1) {
            $control{$data1[0]}{"val"} = 0;
            $control{$data1[0]}{"Gene"} = $data[0];
        }
        if (scalar @data1 ==4) {
            $control{$data1[0]}{"val"} = $data1[1];
            $control{$data1[0]}{"Gene"} = $data[0];
            $control{$data1[0]}{"chr"} = $data1[2];
            $control{$data1[0]}{"exons"} = $data1[3];
        }
    }
    close($fh);
    
    open(my $fh1,"<$infile") or die $!;
    while (my $line = <$fh1>) {
        chomp($line);
        my @data1 = split("\t",$line);
        my @data = split(/\|/,$data1[0]);
        #print STDERR scalar @data."\n";
        if (scalar @data1 ==1) {
            $file_data{$data1[0]}{"val"} = 0;
            $file_data{$data1[0]}{"Gene"} = $data[0];
        }
        if (scalar @data1 ==2) {
            $file_data{$data1[0]}{"val"} = $data1[1];
            $file_data{$data1[0]}{"Gene"} = $data[0];
        }
    }
    close($fh1);
    
    open(my $fho,">$outfile") or die $!;
    
    foreach my $key (keys %control){
        my @key_split = split(/\|/,$key);
        print STDERR "$key\n" unless($key_split[1]);
        if (exists $file_data{$key}) {
            my $tmp = $file_data{$key}{"val"}/$control{$key}{"val"};
            my $pct = sprintf("%.3f",(1-$tmp)*100);
            print $fho $control{$key}{"chr"}."\t".$control{$key}{"Gene"}."\t$key_split[1]\t".$control{$key}{"val"}."\t".$control{$key}{"exons"}."\t$pct%\n";
        }else {
            print $fho $control{$key}{"chr"}."\t".$control{$key}{"Gene"}."\t$key_split[1]\t".$control{$key}{"val"}."\t".$control{$key}{"exons"}."\t100.000%\n";
        }  
    }
    close($fho);
    my $rm_file = $infile;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    
    return $outfile;
}



sub get_pct_perExon{
    
    my ($infile,$control) =@_;
    
    $control=~s/\.txt$/\_perExon\.txt/;
    $infile=~s/\.txt$/\_perExon\.txt/;
    
    unless (-e $control) {print STDERR "$control does not exist\n"; exit;}
    unless (-e $infile) {print STDERR "$infile does not exist\n"; exit;}
    
    my $outfile_pex = $infile;
    $outfile_pex=~s/\_transcriptSIZE\_perExon\.txt/\_perExon\_pct\.txt/;
    
    my %control;
    my %file_data;
    
    open(my $fh,"<$control") or die $!;
    
    while (my $line = <$fh>) {
        chomp($line);
        my @data1 = split("\t",$line);
        #print STDERR $data1[0]."\n";
        if (scalar @data1 ==1) {
            $control{$data1[0]}{"val"} = 0;
        }
        if (scalar @data1 ==5) {
            $control{$data1[0]}{"val"} = $data1[1];
            $control{$data1[0]}{"chr"} = $data1[2];
            $control{$data1[0]}{"start"} = $data1[3];
            $control{$data1[0]}{"stop"} = $data1[4];
        }
    }
    close($fh);
    
    open(my $fh1,"<$infile") or die $!;
    while (my $line = <$fh1>) {
        chomp($line);
        my @data1 = split("\t",$line);
        if (scalar @data1 ==1) {
            $file_data{$data1[0]}{"val"} = 0;
        }
        if (scalar @data1 ==2) {
            $file_data{$data1[0]}{"val"} = $data1[1];
        }
    }
    close($fh1);
    
    open(my $fho,">$outfile_pex") or die $!;
    
    foreach my $key (keys %control){
        my @key_split = split(/\|/,$key);
        my @NM_details = split("_",$key_split[1]);
        if (exists $file_data{$key}) {
            my $tmp = $file_data{$key}{"val"}/$control{$key}{"val"};
            my $pct = sprintf("%.3f",(1-$tmp)*100);
            
            my @regions = `grep -w \"$key\" $combined`;
            my $lowcov_coords='';
            foreach my $lowcov_region (@regions){
                my @tmp = split ("\t",$lowcov_region);
                $lowcov_coords.="$tmp[1]-$tmp[2];";
            }
            $lowcov_coords=~s/;$//;
            if (scalar(@NM_details)==4) {
                print $fho $control{$key}{"chr"}."\t$key_split[0]\t$NM_details[0]_$NM_details[1]\t$NM_details[2]_$NM_details[3]\t".$control{$key}{"start"}."\t".$control{$key}{"stop"}."\t$pct%\t$lowcov_coords\n";
            }else{
                print $fho $control{$key}{"chr"}."\t$key_split[0]\t$NM_details[0]\t$NM_details[1]_$NM_details[2]\t".$control{$key}{"start"}."\t".$control{$key}{"stop"}."\t$pct%\t$lowcov_coords\n";
            }
        }else {
            if (scalar(@NM_details)==4) {
                print $fho $control{$key}{"chr"}."\t$key_split[0]\t$NM_details[0]_$NM_details[1]\t$NM_details[2]_$NM_details[3]\t".$control{$key}{"start"}."\t".$control{$key}{"stop"}."\t100.000%\t.\n";
            }else{
                print $fho $control{$key}{"chr"}."\t$key_split[0]\t$NM_details[0]\t$NM_details[1]_$NM_details[2]\t".$control{$key}{"start"}."\t".$control{$key}{"stop"}."\t100.000%\t.\n";
            }
            
        }   
    }
    close($fho);
    my $rm_file = $infile;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    
    return $outfile_pex;
}


####### Added Feb 5th  #########
sub HGMDcheck {
    my ($pct_file,@genes_lest) = @_;
    my %genes =();
    my $outfile=$pct_file;
    $outfile=~s/-pct\.txt/-pctCov\.txt/;
    foreach my $gene (@genes_lest) {
        chomp($gene);
        $gene=~s/\(/\\\(/;
        $gene=~s/\)/\\\)/;
        my @grep_gene = `grep -w -P \"\t$gene\t\" $pct_file`;
        if ((scalar @grep_gene) == 0) {
            print STDERR "$gene is not present in the provided database\n";
        }
        
        foreach my $gene_trans (@grep_gene){
            chomp($gene_trans);
            my @tmp = split("\t",$gene_trans);
            my $is_HGMD = `grep -w -P \"\t$tmp[2]\$\" $HGMD_db`;
            if ($is_HGMD) {
                $genes{$tmp[1]}{$tmp[2]}{'is_HGMD'} = "true";
            }else{
                $genes{$tmp[1]}{$tmp[2]}{'is_HGMD'} = "false";
            }
            $genes{$tmp[1]}{$tmp[2]}{'line'} = $gene_trans;
            $genes{$tmp[1]}{$tmp[2]}{'cov'} = $tmp[5];
            $genes{$tmp[1]}{$tmp[2]}{'cov'}=~ s/%$//;
        }
    }
    
    open(my $fho," > $outfile") or die $!;
    
    #my $total_genes = keys %genes;
    #print STDERR "$total_genes\n";
    foreach my $gene (sort keys %genes){
        my $written = 0;
        foreach my $transcript (sort keys %{$genes{$gene}}){
            if ($genes{$gene}{$transcript}{'is_HGMD'} eq "true") {
                print $fho "$gene\t$transcript\t$genes{$gene}{$transcript}{'cov'}\tHGMD\n";
                $written=1;
            }
            
        }
        #print STDERR "$gene\t$no_trans\t$written\n";
        if ($written == 0) {
            my $print = "$gene\t";
            my $avg_cov = 0;
            my $no_trans = 0;
            foreach my $transcript (keys %{$genes{$gene}}){
                $no_trans++;
                $avg_cov+=$genes{$gene}{$transcript}{'cov'};
                $print .= "$transcript($genes{$gene}{$transcript}{'cov'});";
            }
            $print=~ s/;$//;
            $avg_cov = sprintf("%0.2f",($avg_cov/$no_trans));
            print $fho "$print\t$avg_cov\n";
        }
    }
    
    close($fho);
}
####### Added Feb 5th  #########

####### Added Feb 10th  #########
sub averageGene {
    my ($pct_file,@genes_lest) = @_;
    my %genes =();
    my $outfile=$pct_file;
    $outfile=~s/-pct\.txt/-pctCov\.txt/;
    foreach my $gene (@genes_lest) {
        chomp($gene);
        $gene=~s/\(/\\\(/;
        $gene=~s/\)/\\\)/;
        my @grep_gene = `grep -w -P \"\t$gene\t\" $pct_file`;
        if ((scalar @grep_gene) == 0) {
            print STDERR "$gene is not present in the provided database\n";
        }
        
        foreach my $gene_trans (@grep_gene){
            chomp($gene_trans);
            my @tmp = split("\t",$gene_trans);
            $genes{$tmp[1]}{$tmp[2]}{'line'} = $gene_trans;
            $genes{$tmp[1]}{$tmp[2]}{'cov'} = $tmp[5];
            $genes{$tmp[1]}{$tmp[2]}{'cov'}=~ s/%$//;
        }
    }
    
    open(my $fho," > $outfile") or die $!;
    
    #my $total_genes = keys %genes;
    #print STDERR "$total_genes\n";
    foreach my $gene (sort keys %genes){
        my $print = "$gene\t";
        my $avg_cov = 0;
        my $no_trans = 0;
        foreach my $transcript (keys %{$genes{$gene}}){
            $no_trans++;
            $avg_cov+=$genes{$gene}{$transcript}{'cov'};
            $print .= "$transcript($genes{$gene}{$transcript}{'cov'});";
        }
        $print=~ s/;$//;
        $avg_cov = sprintf("%0.2f",($avg_cov/$no_trans));
        print $fho "$print\t$avg_cov\n";
    }
    
    close($fho);
   
}
####### Added Feb 10th  #########

##############################################################################################################################

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d", $t->year + 1900, $t->mon + 1, $t->mday, $t->hour, $t->min, $t->sec );
}