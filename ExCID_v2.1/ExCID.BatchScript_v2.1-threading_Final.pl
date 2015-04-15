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
#use Time::localtime;
use Fcntl qw(:flock);
use threads;
use threads::shared;
use Thread::Queue;


unless (scalar @ARGV > 0){USAGE(); exit; }


### GLOBAL VARS ###
# Constant that hold maximum amount of threads to start
use constant MAX_THREADS => 2; 
my %opt;
my $script_dir;
my $script_dir_tmp = $0;
$script_dir_tmp =~m/^.+\//;
$script_dir=$&;
#my $covfasta_gen_script = "java -Xmx8000M -jar $script_dir/bin/CovFasta_Generator.jar ";
my $capcovfasta_gen_script = "java -Xmx8000M -jar $script_dir/bin/CapStatsV2.6.jar ";
my $WGcovfasta_gen_script = "java -Xmx8000M -jar $script_dir/bin/WGSStats_v1.1.jar ";

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
my @bam_files = ();
my @cov_files = ();
my @targeted_database : shared = ();
my $dupRemove = '';

@tmp = `grep "AnnotationDir=" $script_dir/Config.txt`;
$tmp[0]=~s/AnnotationDir=//;
chomp($tmp[0]);
my $annotation_dir = $tmp[0];

my $output_dir = getcwd();
my $viz = '';
my $qc ='';
my $qual ='';
my $batch ='';
my $checkHGMD = '';
my $wholeG ='';
my $bam_list;
my $cov_list;
my @header_cols;
my @low_cov_bedfiles  : shared = ();
my @All_bases_histfiles  : shared = ();
my $GenesCheck_list = ();
my @rm_files  : shared = ();
my $warn = 1;
my $index;
my $is_index_present = -1;
my $annotated_bed_File;
my $LOCK_EXCLUSIVE= 2;
my $UNLOCK= 8;

# A new empty queue
my $q = Thread::Queue->new();
my $q_cov = Thread::Queue->new();

### OPTIONS ###

GetOptions('bam=s' => \@bam_files, 'cov=s' => \@cov_files, 'm:i' => \$opt{m}, 'i:s' => \$opt{i}, 'd' => \$dupRemove, 'qc' => \$qc, 'wig' => \$viz, 'bamList:s' => \$opt{bl}, 'covList:s' => \$opt{cl}, 'p:i' => \$opt{p}, 'qual' => \$qual, 'mq:i' => \$opt{mq}, 'bq:i' => \$opt{bq}, 'batch' => \$batch, 'o:s' => \$opt{o}, 'prefix:s' => \$opt{pre}, 'geneList:s' => \$GenesCheck_list, 'n:i' => \$opt{n}, 'w' => \$wholeG, 'db:s' => \@targeted_database, 'annotatedBed:s' => \$opt{ab}, 'checkHGMD' => \$checkHGMD) || &USAGE;
my $min = $opt{m} || 20; #Default minimum coverage threshold is 20X
my $bed = $opt{i} || "null";
my $pct_samples = $opt{p} || 100;
my $map_q = $opt{mq} || 20;
my $base_q = $opt{bq} || 20;
my $prefix = $opt{pre} || "null";
my $no_of_threads = $opt{n} || 4;

my @bam_files_tmp=();
my @cov_files_tmp=();
if ($opt{bl}) {
    $bam_list = $opt{bl};
    open(my $fh,"< $bam_list") or die $!;
    while (my $bam_file_tm = <$fh>) {
        chomp($bam_file_tm);
        unless(-e $bam_file_tm) {print STDERR "$bam_file_tm from bamlist is not valid. ExCID will skip this file\n";}
        push @bam_files_tmp,$bam_file_tm;
    }
    close($fh);
}

if ($opt{cl}) {
    $cov_list = $opt{cl};
    open(my $fh,"< $cov_list") or die $!;
    while (my $cov_file_tm = <$fh>) {
        chomp($cov_file_tm);
        unless(-e $cov_file_tm) {print STDERR "$cov_file_tm from bamlist is not valid. ExCID will skip this file\n";}
        push @cov_files_tmp,$cov_file_tm;
    }
    close($fh);
}

push @bam_files,@bam_files_tmp;
push @cov_files,@cov_files_tmp;

unless($opt{ab}){
    unless($opt{i}){print STDERR "Please input target BED file. Exiting!!\n"; exit;}

    unless(-e $bed){print STDERR "Please input a valid target BED file. Exiting!!\n"; exit;}
}

unless (@bam_files || @cov_files) {
    print STDERR "No bam files or cov.fasta provided. Exiting!!\n";
    exit;
}

if ($qual) {
    if (!$opt{mq} || !$opt{bq}) {
        print STDERR "No Mapping Quality and Base Quality specified\nUsing 20 as Mapping quality and Base quality threshold.\n";
    }   
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


if (@bam_files) {
    chomp(@bam_files);
    #Enqueue the files
    $q->enqueue(map($_, @bam_files));
    
    # Start the threads and wait for them to finish
    my @threads = map threads->create( \&thread, $q ), 1 .. $no_of_threads;
    $_->join for @threads;
    print STDERR "All Bam Files ran\n" if (scalar(@threads) == scalar(@bam_files));
}
    
if (@cov_files) {
    chomp(@cov_files);

    #Enqueue the files
    $q_cov->enqueue(map($_, @cov_files));
    
    # Start the threads and wait for them to finish
    my @threads = map threads->create( \&thread_covfasta, $q_cov ), 1 .. $no_of_threads;
    $_->join for @threads;
    print STDERR "All Cov Files ran\n" if (scalar(@threads) == scalar(@cov_files));
}


if ($batch) {
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
    
    if ($viz) {
        print STDERR "Generating Visualization track files for Batch analyses...\n\n";
        $sample_list = join(' ', @All_bases_histfiles);
        
        system("mkdir -p $output_dir/split_files/");
        system("rm -rf $output_dir/split_files/*.txt");
        
        foreach my $sample_file (@All_bases_histfiles) {
            system("awk \'{print >> \"$output_dir/split_files/\"\$1\".txt\"}\' $sample_file");
        }
        my $path = "$output_dir/split_files/" ;
        my @files=();
        if (-d $path) {		#Checks if the given path is a directory or not
            @files = <$path/*>;	#Globs all the files in the diven directory
        }else {die "Unable to open $path: $!";}
        
        foreach my $chrom_file (@files) {
            system("$script_dir/bin/bedtools sort -i $chrom_file > $chrom_file-tmp ; mv $chrom_file-tmp $chrom_file");
            system("$script_dir/bin/bedtools groupby -i $chrom_file -g 1,2,3 -c 4 -o collapse >> $output_dir/batch_cov_info.txt");
        }
        
        my %regions;
        open(my $fh,"< $output_dir/batch_cov_info.txt") or die $!;
        open(my $out_avg_hist, "> $output_dir/batch_cov_info.hist") or die $!;
        open(my $out_loss_hist, "> $output_dir/expected_loss.hist") or die $!;
        while (my $line = <$fh>) {
            chomp($line);
            my @data = split(",",$line);
            my @base_coverages;
            for my $i (3 .. scalar(@data)-1){
                push (@base_coverages,int($data[$i]));
            }
            my $average = average(\@base_coverages);
            my $expected_loss = expected_loss(\@base_coverages);
            
            
            print $out_avg_hist "$data[0]\t$data[1]\t$average\n";
            print $out_loss_hist "$data[0]\t$data[1]\t$expected_loss\n";
        }
        close($fh);    
        close($out_avg_hist);
        close($out_loss_hist);
        
        hist_to_WIG ("$output_dir/batch_cov_info.hist");
        hist_to_WIG ("$output_dir/expected_loss.hist");
        
        push @rm_files,"$output_dir/batch_cov_info.hist";
        push @rm_files,"$output_dir/expected_loss.hist";
    }
    
}

push @rm_files, @All_bases_histfiles;

foreach my $waste (@rm_files){
    unlink $waste or warn "Could not unlink $waste: $!";
}
print STDERR "Done\n";


########## SUBROUTINES ###########
sub USAGE {
    print "\nUSAGE: $0 -bam <BAM file> -m <min threshold> -i <Bed file>\n";
    print "  -bam:  BAM file. This option can be used multiple times for more than 1 bam. (Provide full paths to the BAM files)\n";
    print "  -bamList: A text file with 1 Bam file per line.\n";
    print "  -cov:  cov.fasta file. This option can be used multiple times for more than 1 cov.fasta\n";
    print "  -covList: A text file with 1 cov.fasta file per line.\n";
    print "  -m:  Minimal coverage threshold; Default 20\n";
    print "  -i:  Target Bed file. Default WGL Exome3.0\n";
    print "  -d: To discard duplicate reads from coverage calculation. (Duplicate marked Bam required)\n";
    print "  -qc: Provides Gene % coverage and per Exon Coverage for the provided Genelist. Default GeneTest Genes and RefSeq Genes.\n";
    print "  -geneList: List of HGNC approved gene symbols. This options has to be used with -qc . \n";
    print "  -db: A database for gene|transcript. This options has to be used with -qc. (Please check README file for format.) \n";
    print "  -checkHGMD: Check if a Transcript Belongs to HGMD list. If yes then report only that transcript for the gene.\n";
    print "  -batch: Analyses the input files as a Batch. Provides information on bases covered below the coverage threshold in all the files. The Ðp parameter can be used to reduce the stringency.\n";
    print "  -p: % of samples having the base below coverage threshold. Always accompanied by -batch option.\n";
    print "  -wig: Generates Visualization track files. If used with batch option, it will generate the ÒAverage coverage per baseÓ and ÒLoss of coverageÓ tracks.\n";
    print "  -prefix:  Uses this as a prefix for all output files.\n";
    print "  -o: Directory to put all the output files.\n";
    print "  -qual: Provides Quality aware Coverage information where MQ of a read and BQ are as input options.(Default is 20 for both)\n";
    print "  -mq: Sets minimum Mapping Quality of the read to be considered for coverage information. \n";
    print "  -bq: Sets minimum Base Quality to be considered for coverage information.\n";
    print "  -n: maximum number of threads. Default is 4.\n";
    print "  -w: Write Whole Genome cov.fasta file.\n\n";
    exit;
}

sub thread {
    my ($q) = @_;
    
    while (my $bam_file = $q->dequeue_nb()){
        chomp($bam_file);
        next unless($bam_file);
        unless( -e $bam_file){print STDERR "$bam_file does not exist! Please check the file. ExCID will skip this file.\n"; }
        my $out_name="";
        
        if ($prefix ne "null") {
            $out_name= "$prefix";
        }else{
            print STDERR "No prefix added Using BAM file name.\n";
            $out_name = $bam_file;
            $out_name =~ s{.*/}{}; # remove path
            $out_name=~s/\.bam//;
        }
        
        $out_name=$output_dir."/".$out_name;
        
        if ($dupRemove) {
            if ($qual) {
                print STDERR "Generating cov.fasta file, excluding duplicates from calculation...\n\n";
                if ($wholeG) {
                    #system("$covfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -d -m $map_q -b $base_q -w ");
                    system("$WGcovfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -d -m $map_q -b $base_q ");
                    
                }else{
                    #system("$covfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -d -m $map_q -b $base_q ");
                    system("$capcovfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -d -m $map_q -b $base_q ");
                }   
            }else{
                print STDERR "Generating cov.fasta file, excluding duplicates from calculation...\n\n";
                if ($wholeG) {
                    #system("$covfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -d -w ");
                    system("$WGcovfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -d ");
                }else{
                    #system("$covfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -d");
                    system("$capcovfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -d ");
                }   
            }
        }else {
            if ($qual) {
                print STDERR "Generating cov.fasta file...\n\n";
                if ($wholeG) {
                    #system("$covfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -m $map_q -b $base_q -w ");
                    system("$WGcovfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -m $map_q -b $base_q ");
                }else{
                    #system("$covfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -m $map_q -b $base_q ");
                    system("$capcovfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -m $map_q -b $base_q ");
                }
            }else{
                print STDERR "Generating cov.fasta file...\n\n";
                if ($wholeG) {
                    #system("$covfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file -w ");
                    system("$WGcovfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file ");
                }else{
                    #system("$covfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file");
                    system("$capcovfasta_gen_script -o $out_name -t $annotated_bed_File -i $bam_file");
                }
            }
        }
        
        
        print STDERR "Running ExCID...\n\n";
        if ($wholeG) {
            excid_main_WG("$out_name.wholeGenomeCov.fasta",$min,$annotated_bed_File);
        }else{
            excid_main("$out_name.cov.fasta",$min,$annotated_bed_File);
        }
        
        my $low_cov_file = $out_name."_below20x_REPORT.txt";
        push (@low_cov_bedfiles,$low_cov_file);
        
        if ($qc) {
            system("$greping_script -i $annotated_bed_File -l $low_cov_file -db $genetest_DB");
            #system("$greping_script -i $annotated_bed_File -l $low_cov_file -db $refseq_DB");
        }
        if ($GenesCheck_list && -e $GenesCheck_list) {
            if ($checkHGMD) {
                system("$greping_script -i $annotated_bed_File -l $low_cov_file -list $GenesCheck_list -checkHGMD ");
            }else{
                system("$greping_script -i $annotated_bed_File -l $low_cov_file -list $GenesCheck_list ");
            }
        }
        #print STDERR "@targeted_database\n";
        if (scalar(@targeted_database)  > 0) {
            foreach my $target_db (@targeted_database) {
                if ($checkHGMD) {
                    system("$greping_script -i $annotated_bed_File -l $low_cov_file -db $target_db -checkHGMD ");
                }else{
                    system("$greping_script -i $annotated_bed_File -l $low_cov_file -db $target_db ");
                }
            }
        }
        
        if ($batch) {
            my $all_cov_file = $out_name."_AllSites_REPORT.txt_tmp";
            push (@All_bases_histfiles,$all_cov_file);
        }
        
        if ($viz) {
            print STDERR "Generating Visualization track file...\n\n";
            if($wholeG){
                generate_WIG("$out_name.wholeGenomeCov.fasta",$low_cov_file);
            }else{
                generate_WIG("$out_name.cov.fasta",$low_cov_file);
            }
        }  
    }
}

sub thread_covfasta {
    my ($q_cov) = @_;
    
    while (my $cov_file = $q_cov->dequeue_nb()){
        chomp($cov_file);
        next unless($cov_file);
        unless( -e $cov_file){print STDERR "$cov_file does not exist! Please check the file. ExCID will skip this file.\n"; }
        
        my $out_name="";
        if ($prefix ne "null") {
            $out_name= "$prefix";
        }else{
            print STDERR "No prefix added Using COV file name.\n";
            $out_name = $cov_file;
            $out_name =~ s{.*/}{}; # remove path
            if ($wholeG) {
                $out_name=~s/\.wholeGenomeCov\.fasta//;
            }else{
                $out_name=~s/\.cov\.fasta//;
            }
        }
        
        $out_name=$output_dir."/".$out_name;
        
        system("ln -s $cov_file $output_dir/");
        print STDERR "Running ExCID...\n\n";
        if ($wholeG) {
            excid_main_WG("$out_name.wholeGenomeCov.fasta",$min,$annotated_bed_File);
        }else{
            excid_main("$out_name.cov.fasta",$min,$annotated_bed_File);
        }
        my $low_cov_file = $out_name."_below20x_REPORT.txt";
        push (@low_cov_bedfiles,$low_cov_file);
        
        if ($qc) {
            system("$greping_script -i $annotated_bed_File -l $low_cov_file -db $genetest_DB");
            #system("$greping_script -i $annotated_bed_File -l $low_cov_file -db $refseq_DB");
        }
        
        if ($GenesCheck_list && -e $GenesCheck_list) {
            if ($checkHGMD) {
                system("$greping_script -i $annotated_bed_File -l $low_cov_file -list $GenesCheck_list -checkHGMD");
            }else{
                system("$greping_script -i $annotated_bed_File -l $low_cov_file -list $GenesCheck_list");
            }
        }
        
        if (scalar(@targeted_database)  > 0) {
            foreach my $target_db (@targeted_database) {
                if ($checkHGMD) {
                    system("$greping_script -i $annotated_bed_File -l $low_cov_file -db $target_db -checkHGMD");
                }else{
                    system("$greping_script -i $annotated_bed_File -l $low_cov_file -db $target_db");
                }
            }
        }
        
        if ($batch) {
            my $all_cov_file = $out_name."_AllSites_REPORT.txt_tmp";
            push (@All_bases_histfiles,$all_cov_file);
        }
        
        if ($viz) {
            print STDERR "Generating WIG file...\n\n";
            if($wholeG){
                generate_WIG("$out_name.wholeGenomeCov.fasta",$low_cov_file);
            }else{
                generate_WIG("$out_name.cov.fasta",$low_cov_file);
            }
            
        }   
    }
}


sub excid_main_WG{ 
    my ($cov_file, $coverage, $BEDfile) = @_ ;
    my %coord;
    my @bases_covered = ();
    my $total_bases=0;
    my $total_coverage=0;
    my %count;
    my %count_median;
    my $out_name = $cov_file;
    $out_name =~ s{.*/}{}; # remove path
    $out_name=~s/\.wholeGenomeCov\.fasta//;
    
    my $low_bases_out = $output_dir."/".$out_name."_below20x_REPORT.txt";
    my $All_bases_out = $output_dir."/".$out_name."_AllSites_REPORT.txt";
    
    unless(-e $cov_file) {print STDERR "$cov_file does not exist.\n Omitting it from the analyses...";}
        
    my $bed_name = basename($BEDfile);
    $bed_name=~s/\.bed$//;
    
    open(my $FILE,"< $cov_file") or die $!;
    open(my $lcov,"> $low_bases_out") or die "Can't open $low_bases_out:$!\n";
    open(my $hcov,"> $All_bases_out") or die "Can't open $All_bases_out:$!\n";    
    
    my $index = pack 'd', 0;
    $index .= pack 'd', tell $FILE while <$FILE>;
    
    my $size = length(pack("d",0));
    my $i = 0;
    my $line = unpack( 'd', substr($index, $i*$size, $size) );
    my $line_file;
    #print STDERR  "$line \n";
    seek($FILE, $line, 0);
    chomp($line_file = <$FILE> );
    
    while(! eof $FILE) {
        #print "$line_file \n";
        my @tmp = split(" ",$line_file);
        my $chrom = $tmp[0];
        $chrom=~ s/>//;
        $chrom =~ s/[Cc]hr//; #if the number has a 'chr' prefix, remove it.
        my $chrom_size = $tmp[2];
        my $chrom_lines = ceil($chrom_size/100);
        #print "$line_file\t$i\t".($chrom_lines+$i)."\n";
        my @grep_chrom_regions = `grep -w \"^$chrom\" $BEDfile`;
        #print "grep -w \"$chrom\" $BEDfile\n";
        #print STDERR @grep_chrom_regions;
        
        foreach my $chrom_regions (@grep_chrom_regions) {
            chomp($chrom_regions);
            my @chrom_regions_tmp = split("\t",$chrom_regions);
            my $long_coverage ="";
            my @data;
            my $chrom = $chrom_regions_tmp[0];
            my $start = $chrom_regions_tmp[1];
            my $end = $chrom_regions_tmp[2];
            
            
            my $j = $i + ceil($chrom_regions_tmp[1]/100);
            my $j_left = $chrom_regions_tmp[1] % 100;
            my $l = $i + ceil($chrom_regions_tmp[2]/100);
            my $l_left = $chrom_regions_tmp[2] % 100;
            
            $l+=1 if ($l_left == 0);
            $j+=1 if ($j_left == 0);
            
            if ($j == $l) {
                $line = unpack( 'd', substr($index, $j*$size, $size) );
                seek($FILE, $line, 0);
                chomp( $line_file = <$FILE> );
                my @jtmp = split(" ",$line_file);
                
                for(my $k = $j_left ; $k <= $l_left; $k++){
                    $long_coverage.="$jtmp[$k] ";
                    #if ($jtmp[$k] < $coverage) {
                    #    print $lcov "$chrom\t$start\t".($start+1)."\t$jtmp[$k]\n";
                    #}
                    #print $hcov "$chrom\t$start\t".($start+1)."\t$jtmp[$k]\n";
                    #$start++;
                }
                
            }else{
                
                my $length = $chrom_regions_tmp[2] - $chrom_regions_tmp[1];
                my $is_big = 0;
                my $counter=0;
                $is_big = 1 if($length > 10000);
                
                for(my $p = $j ; $p <= $l; $p++){
                    $line = unpack( 'd', substr($index, $p*$size, $size) );
                    seek($FILE, $line, 0);
                    chomp( my $line_file = <$FILE> );
                    my @jtmp = split(" ",$line_file);
                    my $k = 0;
                    $k = $j_left if ($p == $j);
                    my $k_left = scalar(@jtmp) ;
                    $k_left = ($l_left+1) if ($p == $l);
                    
                    for( ; $k < $k_left ; $k++){
                        $long_coverage.="$jtmp[$k] ";
                        $counter++;
                        #if ($jtmp[$k] < $coverage) {
                        #    print $lcov "$chrom\t$start\t".($start+1)."\t$jtmp[$k]\n";
                        #}
                        #print $hcov "$chrom\t$start\t".($start+1)."\t$jtmp[$k]\n";
                        #$start++;
                    }
                }
                
            }
            $long_coverage.="\n";
            chomp($long_coverage);
            my @split = split(/ /,$long_coverage);
            my $p=0;
            foreach my $base ($start .. $end){
                if ($split[$p] < $coverage) {
                    print $lcov "$chrom\t$base\t".($base+1)."\t$split[$p]\n";
                }
                print $hcov "$chrom\t$base\t".($base+1)."\t$split[$p]\n";
                $p++;
            }
        }
        
        $i += $chrom_lines+1;
        
        $line = unpack( 'd', substr($index, $i*$size, $size) );
        seek($FILE, $line, 0);
        chomp( $line_file = <$FILE> ) if(! eof $FILE);
    };
    
    close($FILE);
    close($lcov);
    close($hcov);
    
    hist_to_bed($All_bases_out);
    hist_to_bed($low_bases_out);
    
}

sub excid_main{ 

    my ($cov_file, $coverage, $BEDfile) = @_ ;
    my %coord;
    my @bases_covered = ();
    my @bases_covered_bed = ();
    my $total_bases=0;
    my $total_coverage=0;
    my %count;
    my %count_median;
    my $out_name = $cov_file;
    $out_name =~ s{.*/}{}; # remove path
    $out_name=~s/\.cov\.fasta//;
    
    my $low_bases_out = $output_dir."/".$out_name."_below20x_REPORT.txt";
    my $All_bases_out = $output_dir."/".$out_name."_AllSites_REPORT.txt";
    
    unless(-e $cov_file) {print STDERR "$cov_file does not exist.\n Omitting it from the analyses...";}
        
    my $bed_name = basename($BEDfile);
    $bed_name=~s/\.bed$//;
    
    open(my $lcov,"> $low_bases_out") or die "Can't open $low_bases_out:$!\n";
    open(my $hcov,"> $All_bases_out") or die "Can't open $All_bases_out:$!\n";    
    open(my $fh,"<$cov_file") or die $!;
    
    my $start;
    my $end;
    my $chrom;
    my $long_coverage;
    while (my $line = <$fh>) {
        chomp($line);
        if ($line=~m/^>/) {
             my @tmp = split(" ",$line);
             $chrom = $tmp[0];
             $chrom=~s/^>//;
             $start = $tmp[1];
             $end = $tmp[2];
             $long_coverage="";
             #print STDERR "$chrom $start $end\n";
             next;
        }
        $long_coverage = $line;
        my @split = split(/ /,$long_coverage);
           
        while (scalar @split != ($end-$start+1)){
            $line = <$fh>;
            chomp($line);
            $long_coverage.= "$line";
            @split = split(/ /,$long_coverage);
        }
        
        my $i=0;
        foreach my $base ($start .. $end){
            my $key = "$chrom"."_$base";
            $coord{"$key"} = $split[$i];
            $i++;
        }
    }
    close($fh);
    
    open(my $fhb,"< $BEDfile") or die $!;
    while (my $line = <$fhb>) {
        chomp($line);
        my @tmp = split("\t",$line);
        $chrom = $tmp[0];
        $start = $tmp[1];
        $end = $tmp[2];
        my @data;
        foreach my $base ($start .. $end){
            my $key = "$chrom"."_$base";
            $total_bases++;
            unless (exists $coord{$key}){
                print STDERR "$chrom\t$base\t".($base+1)."\t0\n";
                $coord{$key} = 0;
                #next;
            }
            
            if ($coord{$key} < $coverage) {
                print $lcov "$chrom\t$base\t".($base+1)."\t$coord{$key}\n";
            }
            print $hcov "$chrom\t$base\t".($base+1)."\t$coord{$key}\n";
        }	
    }
    
    close($fhb);
    close($lcov);
    close($hcov);
    hist_to_bed($All_bases_out);
    hist_to_bed($low_bases_out);
}

sub hist_to_bed {
    my $infile = $_[0]; #histogram
    my $start;
    my $end;
    my $count;
    my $chr_start;
    my @coverage=();
    
    my $outfile = $infile."_tmp";
    
    open(my $INFILE, "< $infile ") || die "Can't open $infile $!\n";
    open(my $OUTFILE, "> $outfile ") || die "Can't open $outfile $!\n";
    
    open (HEADER, "< $index/HEADER.txt") or die "Header file missing in $index $!\n";
    my $header = <HEADER>;
    @header_cols = split(" ",$header);
    my $out_header = join("\t",@header_cols[0..2])."\tLength\tCoverage\t".join("\t",@header_cols[3..(scalar(@header_cols)-1)]);
    #print $OUTFILE '##generatedBy:', "$0", "\n";
    #print $OUTFILE "##fileDate=".datestamp()."\n";
    #print $OUTFILE "##$out_header\n";
    
    
    while (<$INFILE>) {
        chomp $_;
        my ($chr, $coord, $coord_end, $cov) = split(/\t/, $_);
        if (! $count) {
            $chr_start = $chr;
            $start = $coord;
            $count = $coord + 1;
            push @coverage,$cov;
        } elsif ($count != $coord) {
            $end = $count - 1;
            #my $length = $end - $start; # previous
            my $length = $end - $start+1;
            @coverage = sort { $a <=> $b } @coverage;
            my $median = int(@coverage % 2 ? $coverage[(@coverage-1)/2] : ($coverage[@coverage/2-1]+$coverage[@coverage/2])/2);
            my $total = 0;
            foreach my $v (@coverage) {
                $total += $v;
            }
            my $mean = sprintf("%.0f",$total/@coverage);
            #if ($length > 0) { # previous
            if ($length >= 0) {
                #my $anno_result_out = &ANNO_TARGET($index, $chr_start, $start, $end, $length, $mean);
                #print $OUTFILE "$anno_result_out\n";
                print $OUTFILE "$chr_start\t$start\t$end\t$length\t$mean\n";
            }
            $chr_start = $chr;
            $start = $coord;
            $count = $coord + 1;
            @coverage = ();
            push @coverage,$cov;
        } elsif ($count == $coord) {
            $chr_start = $chr;
            $count++;
            push @coverage,$cov;
        } else {
            print "\n\nFail; lapse in logic!!\n\n";
            exit;
        }
    }
    $end = $count - 1;
    
    @coverage = sort { $a <=> $b } @coverage;
    #my $median = int(@coverage % 2 ? $coverage[(@coverage-1)/2] : ($coverage[@coverage/2-1]+$coverage[@coverage/2])/2);
    my $total = 0;
    foreach my $v (@coverage) {
            $total += $v;
    }
    my $mean = sprintf("%.0f",$total/@coverage);
    
    my $length = scalar(@coverage);
    #my $anno_result_out = &ANNO_TARGET($index, $chr_start, $start, $end, $length, $mean);
    #print $OUTFILE "$anno_result_out\n";
    print $OUTFILE "$chr_start\t$start\t$end\t$length\t$mean\n";
    
    close($INFILE);
    close($OUTFILE);
    
    my $cmd = "\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5";
    #for(my $i =3; $i <= scalar(@header_cols) ; $i++){
    for(my $i =3; $i < scalar(@header_cols) ; $i++){  ##check Feb5th
        my $tmp = $i+3+3;
        $cmd .= "\"\\t\"\$$tmp";
    }
    
    system ("$script_dir/bin/bedtools intersect -header -wao -a $outfile -b $annotated_bed_File | awk -F\$\'\t\' '{print $cmd }'  | $script_dir/bin/bedtools sort -i  | sed 's/\t-1\t/\t.\t/'   | uniq    > $outfile-tmp ");
    #system ("mv $outfile-tmp $outfile ");
    
    my $tmp_header = "##generatedBy:'".$0."\n";
    $tmp_header .=  "##fileDate=".datestamp()."\n";
    $tmp_header .=  "##TargetsBedFile=$bed\n";
    $tmp_header .=  "##$out_header";
    
    system ("echo \"$tmp_header\" > $outfile");
    system ("cat $outfile-tmp >> $outfile");
    
    #my $rm_file = "$outfile-tmp" ;
    #unlink $rm_file or warn "Could not unlink $rm_file: $!";
    push @rm_files, "$outfile-tmp" ;
    
    
    if ($batch && $infile=~m/AllSites_REPORT/) {
        system("cp $outfile $outfile-tmp");
        system("mv $infile $outfile");
        system("mv $outfile-tmp $infile");
    }else{
        system("mv $outfile $infile");
    }
    
    
}


sub ANNO_TARGET {
    my $trigger = 11000; #The biggest exon I've seen is ~10k, no use looking beyond that.
    my $index = shift;
    my $lowchr = shift;
    my $lowstart = shift;
    my $lowstop = shift;
    my $length = shift;
    my $cov = shift;
    my $annoresult;

    my $index_path = "$index/";
    my $index_file = $index_path . $lowchr . ".targetindex.txt"; 
    open (INDEX, "<$index_file") || die "Can't open $index_file: $!\n";
    INDEXLOOP: while (<INDEX>) {
	    chomp $_;
	    my $line = $_;
	    my @columns = split(/\t/, $line);
	    if ($warn == 1 && scalar @columns != scalar @header_cols ) {
		print STDERR "The number of columns in Header do not match the index files. Please check the Header.txt file.\n";
		$warn = 0;
	    }
	    
	    my $start = $columns[1];
	    my $stop = $columns[2];
	    my $start_dif = $lowstart - $start;
	    if ( $start_dif >= $trigger ) {
	        next INDEXLOOP;
	    } elsif ( ($start <= $lowstart) && ($stop >= $lowstop)) {
                if ($cov eq "null") {
                    $annoresult = "$lowchr\t$lowstart\t$lowstop\t$length\t ";
                }else{
                    $annoresult = "$lowchr\t$lowstart\t$lowstop\t$length\t$cov";
                }
		for(my $i = 3; $i<scalar(@columns);$i++){
		    $annoresult .="\t$columns[$i]";
		}
		return $annoresult;
	    }
    }
    if (! $annoresult) {
        if ($cov eq "null") {
            $annoresult = "$lowchr\t$lowstart\t$lowstop\t$length\t ";
        }else{
            $annoresult = "$lowchr\t$lowstart\t$lowstop\t$length\t$cov";
        }
        for(my $i = 5; $i<(scalar(@header_cols)+2);$i++){
            $annoresult .="\t.";
        }
	return $annoresult;
    }
    close INDEX;
}


sub generate_WIG {
    
    my ($infile,$lowcov_in) = @_;
    
    (my $outfile = $infile) =~ s{\.[^.]+$}{}; # removes extension
    $outfile=~ s{\.[^.]+$}{}; # removes extension
    $outfile .= "_UCSCtrack.wig";
    
    # All output files will have this header:
    open(OUTFILE, ">$outfile") || die "Can't open $outfile: $!\n";
    print OUTFILE "track type=wiggle_0 name=Coverage autoScale=off viewLimits=0.0:100.0 graphType=bar visibility=full\n";
    
    
    open(INFILE, "<$infile") || die "Can't open $infile: $!\n";
    while (<INFILE>) {
        chomp $_;
        if ($_ =~ m/^\>/) {
            my ($chr, $start_coord, $stop_coord) = split(/ /, $_);
            $chr =~ s/\>//;
            $chr =~ s/[Cc]hr//; #I need to control how 'chr' looks.
            print OUTFILE "fixedStep  chrom=chr$chr  start=$start_coord  step=1\n";
            next;
        } elsif ($_ =~ m/^$/) {
            #this line is blank; skip it.
            next; 
        }
    
        my @array = split(/ /, $_);
    
        foreach (@array) {
                my $num = $_; chomp $_;
                print OUTFILE "$num\n";
        }
    }
    close OUTFILE;
    close INFILE;

    my $lowcov_out = $lowcov_in;
    
    $lowcov_out =~s/\_REPORT\.txt/\_UCSCtrack\.bed/ ;
    
    open(LOWOUT, ">$lowcov_out") || die "Can't open $lowcov_out: $!\n";
    open(LOWIN, "<$lowcov_in") || die "Can't open $lowcov_in: $!\n";
    print LOWOUT "track name=\'Low Coverage\' description=\'Coverage <".$min."X\' color=255,0,0\n";
    while (<LOWIN>) {
        chomp $_;
        if ($_ =~ m/^\#/) {next;}
        if ($_ =~ m/^[0-9]$/) {next;}
        my @low_array = split(/\t/, $_);
        my $lowchr = $low_array[0];
        my $lowstart = $low_array[1];
        my $lowstop = $low_array[2];
        print LOWOUT "chr$lowchr\t$lowstart\t$lowstop\n";
    }
    close LOWOUT;
    close LOWIN;
}


sub average { 
    @_ == 1 or die ('Usage: $average = average(\@array);'); 
    my ($array_ref) = @_; 
    my $sum; 
    my $count = scalar @$array_ref; 
    foreach (@$array_ref) { $sum += $_; }
    my $ret = sprintf "%.2f", ($sum / $count);
    return $ret; 
}


sub expected_loss { 
    @_ == 1 or die ('Usage: $expected_loss = expected_loss(\@array);'); 
    my ($array_ref) = @_; 
    my $expected_loss; 
    my $count = scalar @$array_ref; 
    foreach (@$array_ref) {
        if ($_ < 20) {
            $expected_loss++;
        }
    }
    return $expected_loss; 
}


sub hist_to_WIG {
    my $infile = $_[0]; #histogram
    my $start;
    my $end;
    my $count;
    my $chr_start;
    my @coverage=();
    my $outfile = $infile;
    $outfile=~s/\.hist/\.wig/;
    open(INFILE, "<$infile") || die "Can't open $infile: $!\n";
    open(OUTFILE, "> $outfile") || die "Can't open $infile: $!\n";
    while (<INFILE>) {
        chomp $_;
        my ($chr, $coord, $cov) = split(/\t/, $_);
        if (! $count) {
            $start = $coord;
            $count = $coord + 1;
            push @coverage,$cov;
        } elsif ($count != $coord) {
            $end = $count - 1;
            print OUTFILE "fixedStep  chrom=chr$chr_start  start=$start  step=1\n";
            my $i = 0;
            foreach my $base ($start .. $end){
                print OUTFILE "$coverage[$i]\n" or die;
                $i++;
            }
            $start = $coord;
            $count = $coord + 1;
            @coverage = ();
            push @coverage,$cov;
        } elsif ($count == $coord) {
            $chr_start = $chr;
            push @coverage,$cov;
            $count++;
        } else {
            print STDERR "\n\nFail; lapse in logic!!\n\n";
            exit;
        }
    }
    $end = $count - 1;
    
    print OUTFILE "fixedStep  chrom=chr$chr_start  start=$start  step=1\n";
    my $i = 0;
    foreach my $base ($start .. $end){
        print OUTFILE "$coverage[$i]\n" or die;
        $i++;
    }

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
