#!/usr/bin/perl -w
#Author: Rashesh Sanghvi
#Purpose: ExCID Report serves two main purposes Ð the identification of regions below 20X with base-pair granularity, and the subsequent gene/transcript/exon annotations of those regions.
#Date: 
####################################
use strict;
use warnings;
use diagnostics;
use File::Basename;
use POSIX;
use Getopt::Long;
use File::Path;
use Fcntl qw(:flock);


unless (scalar @ARGV > 0){USAGE(); exit; }


### GLOBAL VARS ###
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
my $genetest_DB = "$database_dir/genetest-Nov2013.bed";
my $refseq_DB = "$database_dir/RefSEQ-Apr2014.bed";


my $greping_script = "$script_dir/ExCID.grep_gene_list_pct_Final.pl";
my @targeted_database = ();
my $dupRemove = '';

@tmp = `grep "AnnotationDir=" $script_dir/Config.txt`;
$tmp[0]=~s/AnnotationDir=//;
chomp($tmp[0]);
my $annotation_dir = $tmp[0];

my $output_dir = getcwd();
my $qc ='';
my $checkHGMD = '';
my $bam_list;
my @header_cols;
my $report_list;
my @low_cov_bedfiles = ();
my $GenesCheck_list = ();
my @rm_files = ();
my $index;
my $is_index_present = -1;
my $annotated_bed_File;
my $LOCK_EXCLUSIVE= 2;
my $UNLOCK= 8;


### OPTIONS ###

GetOptions('i:s' => \$opt{i}, 'qc' => \$qc, 'ReportList:s' => \$report_list, 'p:i' => \$opt{p}, 'o:s' => \$opt{o}, 'geneList:s' => \$GenesCheck_list, 'db:s' => \@targeted_database, 'annotatedBed:s' => \$opt{ab}, 'checkHGMD' => \$checkHGMD) || &USAGE;
my $min = $opt{m} || 20; #Default minimum coverage threshold is 20X
my $bed = $opt{i} || "null";
my $pct_samples = $opt{p} || 100;


if (-e $report_list) {
    open(my $fh,"< $report_list") or die $!;
    while (my $report_file_tm = <$fh>) {
        chomp($report_file_tm);
        unless(-e $report_file_tm) {print STDERR "$report_file_tm from list is not valid. ExCID will skip this file\n";}
        push @low_cov_bedfiles,$report_file_tm;
    }
    close($fh);
}else{
    print STDERR "Please input Valid file $report_list. Exiting!!\n"; exit;
}

unless($opt{ab}){
    unless($opt{i}){print STDERR "Please input target BED file. Exiting!!\n"; exit;}

    unless(-e $bed){print STDERR "Please input a valid target BED file. Exiting!!\n"; exit;}
}

unless (@low_cov_bedfiles) {
    print STDERR "No _below20x_REPORT.txt files provided. Exiting!!\n";
    exit;
}

if ($opt{o}) {
    $output_dir = $opt{o}; 
}


{
    if($opt{ab}){
        unless(-e $opt{ab}){print STDERR "Please input a valid Annotated BED file. Exiting!!\n"; exit;}
        my $tmp_bed_file_name = $opt{ab};
        $tmp_bed_file_name =~ s{.*/}{}; # remove path
        
        my $output_dir_bed = "$annotation_dir/$tmp_bed_file_name-ouput_dir";
        $index = "$output_dir_bed/index_files";
        $annotated_bed_File = "$output_dir_bed/$tmp_bed_file_name.Annotated";
        $bed = "$output_dir_bed/$tmp_bed_file_name.Annotated";
        system("mkdir -p $index");
        system("cp $opt{ab} $annotated_bed_File");
        #system("awk \'{print > \"$index/\"\$1\".targetindex.txt\"}\' $annotated_bed_File");
        system("awk \'{print \$0 >> (\"$index/\"\$1\".targetindex.txt\") ; close(\"$index/\"\$1\".targetindex.txt\")}\' $annotated_bed_File");
        open(my $fhb,"<$annotated_bed_File") or die $!;
        my $line = <$fhb>;
        chomp($line);
        my @tmp = split("\t",$line);
        my $cols = scalar @tmp;
        close($fhb);
        
        my $header = "Chr\tStart\tStop";
        foreach my $i (3 .. ($cols-1)){
            my $tmp_cnt = $i-6;
            $header .= "\tUsr_col_$tmp_cnt";
        }
        $header .="\n";
        open(my $h, ">$index/HEADER.txt") or die $!;
        
        print $h $header;
        close($h);
        chmod(0775, $output_dir_bed) or print STDERR "Couldn't chmod $output_dir_bed: $!";
        goto INDEXPRESENT;
    }
}

{
    my $tmp_bed_file_name = $bed;
    $tmp_bed_file_name =~ s{.*/}{}; # remove path
    
    my $output_dir_bed = "$annotation_dir/$tmp_bed_file_name-ouput_dir";
    $index = "$output_dir_bed/index_files";
    
    my $SEMAPHORE = $output_dir_bed.".lck";
    my $time = 0;
    my $SEMlock;
    my $random_number = rand(100);
    sleep($random_number);
    do{
        goto Loop if ($time >= 90);
        my $lock_status=check_lock_exists($SEMAPHORE);
        if ($lock_status == 0){
            print "-----Another instance is accessing the $output_dir_bed, sleeping for 10 seconds------\n";
            $time+=10;
            sleep(10);
        }else{
            print STDERR "Lock file $SEMAPHORE\n";
            open ($SEMlock, ">$SEMAPHORE") || die "$!\nproblem opening lock File\n";
            flock $SEMlock, $LOCK_EXCLUSIVE;
            print "----Obtained Lock----".timestamp()."\n";
        }
    }while(!$SEMlock);
    
    $time = 0;
    
    Loop: {
        $time = 0;
        do{
            if (-d $output_dir_bed) {
                if (-d $index) {
                    if (-e "$index/HEADER.txt") {
                        print STDERR "The index files are already present for $tmp_bed_file_name\n";
                        $is_index_present = 1;
                        $annotated_bed_File = "$output_dir_bed/$tmp_bed_file_name.Annotated";
                        flock $SEMlock, $UNLOCK;
                        close($SEMlock);
                        unlink($SEMAPHORE) or die "can't remove lockfile: $SEMAPHORE ($!)".timestamp()."\n";
                        print "----Unlocked----\n";
                        goto INDEXPRESENT;
                    }else{
                        sleep(5);
                        $time+=5;
                        last Loop if ($time >= 60);
                    }
                }else{
                    sleep(5);
                    $time+=5;
                    last Loop if ($time >= 60);
                }
            }else{last Loop ;} 
        }while($is_index_present != 1);
    }
    
    
    system("mkdir -p $output_dir_bed");
    
    opendir(my $Dir, "$output_dir_bed/") or die "Can't open $output_dir_bed ($!)";
    
    print STDERR "Annotating Bed file...\n\n";
    
    system("ln -s $bed $output_dir_bed/");
    
    $bed =~ s{.*/}{}; # remove path
    $bed = "$output_dir_bed/$bed";
    $annotated_bed_File = "$bed.Annotated";
    system ("perl $script_dir/bed_file-annotator_V2_RefSeq-VEGA.pl $bed $database1 $output_dir_bed/" );
    system ("perl $script_dir/bed_file-annotator_V2_CCDS-miRBASE.pl $bed $database2 $output_dir_bed/" );
    system ("perl $script_dir/bed_file-annotator_V2_RefSeq-VEGA.pl $bed $database3 $output_dir_bed/" );
    system ("perl $script_dir/bed_file-annotator_V2_CCDS-miRBASE.pl $bed $database4 $output_dir_bed/" );
    
    my $file_tmp1 = $bed."-".basename($database1);
    my $file_tmp2 = $bed."-".basename($database2);
    my $file_tmp3 = $bed."-".basename($database3);
    my $file_tmp4 = $bed."-".basename($database4);
    
    open(my $fhb,"<$bed") or die $!;
    my $line = <$fhb>;
    chomp($line);
    my @tmp = split("\t",$line);
    my $cols = scalar @tmp;
    close($fhb);    
    
    my $cmd;
    open($fhb,"<$file_tmp1") or die $!;
    $line = <$fhb>;
    chomp($line);
    @tmp = split("\t",$line);
    my $cols1 = scalar @tmp;
    close($fhb);
    
    open($fhb,"<$file_tmp2") or die $!;
    $line = <$fhb>;
    chomp($line);
    @tmp = split("\t",$line);
    my $cols2 = scalar @tmp;
    close($fhb);
    
    open($fhb,"<$file_tmp3") or die $!;
    $line = <$fhb>;
    chomp($line);
    @tmp = split("\t",$line);
    my $cols3 = scalar @tmp;
    close($fhb);
    
    open($fhb,"<$file_tmp4") or die $!;
    $line = <$fhb>;
    chomp($line);
    @tmp = split("\t",$line);
    my $cols4 = scalar @tmp;
    close($fhb);
    
    $cmd="\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$".(1+$cols1+3)."\"\\t\"\$".(1+$cols1+$cols2+1+3)."\"\\t\"\$".(1+$cols1+$cols2+1+$cols3+1+3)."\"\\t\"\$6\"\\t\"\$".(1+$cols1+3+2)."\"\\t\"\$".(1+$cols1+$cols2+1+3+2)."\"\\t\"\$".(1+$cols1+$cols2+1+$cols3+1+3+2)."\"\\t\"\$5\"\\t\"\$".(1+$cols1+3+1)."\"\\t\"\$".(1+$cols1+$cols2+1+3+1)."\"\\t\"\$".(1+$cols1+$cols2+1+$cols3+1+3+1)."";    
    
    for(my $i =1; $i <= ($cols-3); $i++){
        my $tmp = (1+$cols1+$cols2+1+$cols3+1+3+2)+$i;
        $cmd .= "\"\\t\"\$$tmp";
    }
    
    system ("$script_dir/bin/bedtools intersect -wao -a $file_tmp1 -b $file_tmp2 | $script_dir/bin/bedtools intersect -wao -a - -b $file_tmp3 | $script_dir/bin/bedtools intersect -wao -a - -b $file_tmp4  | awk -F\$\'\t\' '{print $cmd }' > $annotated_bed_File-tmp");
    
    system ("perl $script_dir/reformat.pl $annotated_bed_File-tmp > $annotated_bed_File ");
    
    open($fhb,"<$annotated_bed_File") or die $!;
    $line = <$fhb>;
    chomp($line);
    @tmp = split("\t",$line);
    $cols = scalar @tmp;
    close($fhb);
    
    my $header = "Chr\tStart\tStop\tGene\tPrev\tSynonymous\tRefSeq\tCCDS\tVEGA\tmiRNA";
    foreach my $i (10 .. ($cols-1)){
        my $tmp_cnt = $i-6;
        $header .= "\tUsr_col_$tmp_cnt";
    }
    $header .="\n";
    system("mkdir -p $index");
    
    system("awk \'{print \$0 >> (\"$index/\"\$1\".targetindex.txt\") ; close(\"$index/\"\$1\".targetindex.txt\")}\' $annotated_bed_File");
    
    my $rm_file = "$annotated_bed_File-tmp";
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = $file_tmp1 ;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = $file_tmp2 ;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = $file_tmp3 ;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = $file_tmp4 ;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    
    unlink($bed);
    
    open(my $h, ">$index/HEADER.txt") or die $!;
    
    print $h $header;
    close($h);
    closedir($Dir);
    
    flock $SEMlock, $UNLOCK;
    close($SEMlock);
    unlink($SEMAPHORE) or die "can't remove lockfile: $SEMAPHORE ($!)";
    print "----Unlocked----\n";
    
    chmod(0775, $output_dir_bed) or print STDERR "Couldn't chmod $output_dir_bed: $!";
}
INDEXPRESENT:

{
    my $samples_to_check = scalar(@low_cov_bedfiles);
    
    my $sample_numbers_to_check = int(($pct_samples/100)*($samples_to_check));
    my $sample_list = join(' ', @low_cov_bedfiles);
    
    system("$script_dir/bin/bedtools multiinter -i $sample_list > $output_dir/batch_Lowcov.bed");
    
    open(my $out_batch, "< $output_dir/batch_Lowcov.bed") or die $!;
    open(my $out_batch_threshold, "> $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed") or die $!;
    while (my $line = <$out_batch>) {
        chomp ($line);
        my @split = split("\t",$line);
        if ($split[3] >= $sample_numbers_to_check) {
            my $length = $split[2]-$split[1];
            #my $anno_result_out = &ANNO_TARGET($index, $split[0], $split[1], $split[2], $length,"null");
            #print $out_batch_threshold "$anno_result_out\n";
            print $out_batch_threshold "$split[0]\t$split[1]\t$split[2]\t$length\n";
        }
    }
    close($out_batch);
    close($out_batch_threshold);
    
    open (HEADER, "< $index/HEADER.txt") or die "Header file missing in $index $!\n";
    my $header = <HEADER>;
    @header_cols = split(" ",$header);
    my $out_header = join("\t",@header_cols[0..2])."\tLength\t".join("\t",@header_cols[3..(scalar(@header_cols)-1)]);
    close(HEADER);
    
    my $cmd = "\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4";
    for(my $i =3; $i < scalar(@header_cols) ; $i++){
        my $tmp = $i+3+2;
        $cmd .= "\"\\t\"\$$tmp";
    }
    
    system ("$script_dir/bin/bedtools intersect -wao -a $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed -b $annotated_bed_File | awk -F\$\'\t\' '{print $cmd }' > $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed-tmp");
    #system ("mv $outfile-tmp $outfile ");
    
    my $tmp_header = "##generatedBy:'".$0."\n";
    $tmp_header .=  "##fileDate=".datestamp()."\n";
    $tmp_header .=  "##TargetsBedFile=$bed\n";
    $tmp_header .=  "##$out_header";
    
    system ("echo \"$tmp_header\" > $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed");
    system ("cat $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed-tmp >> $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed");
    
    #my $rm_file = "$output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed-tmp" ;
    #unlink $rm_file or warn "Could not unlink $rm_file: $!";
    
    push @rm_files, "$output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed-tmp" ;
    push @rm_files, "$output_dir/batch_Lowcov.bed";
    
    if ($qc) {
        system("$greping_script -i $annotated_bed_File -l $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed -db $genetest_DB");
        #system("$greping_script -i $annotated_bed_File -l $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed -db $refseq_DB");
    }
    if ($GenesCheck_list && -e $GenesCheck_list) {
        if ($checkHGMD) {
            system("$greping_script -i $annotated_bed_File -l $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed -list $GenesCheck_list -checkHGMD");
        }else{
            system("$greping_script -i $annotated_bed_File -l $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed -list $GenesCheck_list ");
        }
    }
      
    if (scalar(@targeted_database)  > 0) {
        foreach my $target_db (@targeted_database) {
            if ($checkHGMD) {
                system("$greping_script -i $annotated_bed_File -l $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed -db $target_db -checkHGMD");
            }else{
                system("$greping_script -i $annotated_bed_File -l $output_dir/batch_$sample_numbers_to_check-samples_Lowcov.bed -db $target_db");
            }
        }
    } 
    
}

foreach my $waste (@rm_files){
    unlink $waste or warn "Could not unlink $waste: $!";
}
print STDERR "Done\n";


########## SUBROUTINES ###########
sub USAGE {
    print "\nUSAGE: $0 -ReportList <List of low cov files> -i <Bed file> -p <% of samples>\n";
    print "  -ReportList: A list of *_below20x_REPORT.txt. One file per line. Please provide full paths.\n";
    print "  -i:  Target Bed file. Default WGL Exome3.0\n";
    print "  -qc: Provides Gene % coverage and per Exon Coverage for the provided Genelist. Default GeneTest Genes and RefSeq Genes.\n";
    print "  -geneList: List of HGNC approved gene symbols. This options has to be used with -qc . \n";
    print "  -db: A database for gene|transcript. This options has to be used with -qc. (Please check README file for format.) \n";
    print "  -checkHGMD: Check if a Transcript Belongs to HGMD list. If yes then report only that transcript for the gene.\n";
    print "  -p: % of samples having the base below coverage threshold. Always accompanied by -batch option.\n";
    print "  -o: Directory to put all the output files.\n";
    exit;
}

##############################################################################################################################

sub timestamp {
  my $t = localtime;
  return $t;
}

sub datestamp {
  my $t = localtime;
  my @tmp = split(" ", $t);
  return "$tmp[4]-$tmp[1]-$tmp[2]";
}

#Sub Routine to Check Whether the Lock File Exists
sub check_lock_exists {
  
   my $lock_file = $_[0];
   my $status;

   if (-e $lock_file) {
        $status = 0;
    } else {
        $status = 1;
    }
  
   return $status;
 }
