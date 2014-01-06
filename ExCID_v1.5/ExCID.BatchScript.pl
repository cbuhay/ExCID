#!/usr/bin/perl -w
## Wraper Script for the ExCID.
use strict;
use diagnostics;
use Getopt::Long;



unless (scalar @ARGV > 0){USAGE(); exit; }


### GLOBAL VARS ###
my %opt;
my $script_dir;
my $script_dir_tmp = $0;
$script_dir_tmp =~m/^.+\//;
$script_dir=$&;
my $index = "$script_dir/index_files";
my $output_dir = "$script_dir/ouput_dir";
my $database1 = "$script_dir/database/whole_UCSC_detail.txt.final-noCHR.bed";
my $database2 = "$script_dir/database/whole_Refseq_UCSC_detail.txt.final-noCHR.bed";
my $script ="$script_dir/ExCID_getlowcov_and_annotate_exome_v1.5.pl";
my $default_index = "$script_dir/VCRome_index_131206/";
my $covfasta_gen_script = "java -Xmx6000M -jar $script_dir/bin/Generate_covfasta.jar";

### OPTIONS ###
my @bam_files = ();
GetOptions('f=s' => \@bam_files, 'm:i' => \$opt{m}, 'i:s' => \$opt{i}) || &USAGE;
my $min = $opt{m} || 20; #meant for sites <20X coverage
my $bed = $opt{i} || "null";
my $file;


if ($bed ne "null") {
    print STDERR "Annotationg Bed file...\n\n ";
    system("mkdir -p $output_dir");
    system("ln -s $bed $output_dir/");
    
    $bed =~ s{.*/}{}; # remove path
    $bed = "$output_dir/$bed";
    my %bedfile;
    
    open(my $fhb,"<$bed") or die $!;
    my $line = <$fhb>;
    chomp($line);
    my @tmp = split("\t",$line);
    my $cols = scalar @tmp;
    close($fhb);
    
    my $cmd="";
    
    for(my $i =0; $i < $cols; $i++){
        my $tmp = $i+1;
        $cmd .= "\$$tmp\"\t\"";
    }
    $cmd .= "\$".($cols+3+1)."\t";
    
    system("$script_dir/bin/bedtools intersect -a $bed -b $database1 -wao | awk -F\$\'\t\' '{print $cmd}' > $bed.Annotated-tmp");
    
    
    
    open($fhb,"<$bed.Annotated-tmp") or die $!;
    $line = <$fhb>;
    chomp($line);
    @tmp = split("\t",$line);
    $cols = scalar @tmp;
    close($fhb);
    
    $cmd="";
    for(my $i =0; $i < $cols; $i++){
        my $tmp = $i+1;
        $cmd .= "\$$tmp\"\t\"";
    }
    $cmd .= "\$".($cols+3+1)."\t";
    
    system("$script_dir/bin/bedtools intersect -a $bed.Annotated-tmp -b $database2 -wao | awk -F\$\'\t\' '{print $cmd}' > $bed.Annotated");
    system("mkdir -p $index");
    
    $file = "$bed.Annotated";
    my %data_annotation;
    open(my $fh, "<$file") or die $!;
    while (my $line = <$fh>) {
        chomp($line);
        my @columns = split("\t",$line);
        my $key = "$columns[0]_$columns[1]_$columns[2]";
        unless (exists $data_annotation{$key}){
            $data_annotation{$key}{'gene'} = "";
            $data_annotation{$key}{'transcript'} = "";
            $data_annotation{$key}{'other'} = "";
            $data_annotation{$key}{'rest'} = join("\t", @columns[3..(scalar(@columns)-3)]);
        }
        
        my $db1 = $columns[scalar(@columns)-2];
        my $db2 = $columns[scalar(@columns)-1];
        my @db1_columns = split(";",$db1);
        my @db2_columns = split(";",$db2);
        
        foreach my $db1_col (@db1_columns){
            if ($db1_col eq ".") {
                last;
            }
            
            my @db1_cols = split(/\|/,$db1_col);
            if (length($data_annotation{$key}{'gene'}) == 0) {
                $data_annotation{$key}{'gene'} = $db1_cols[0];
            }else {
                unless (index($data_annotation{$key}{'other'}, $db1_cols[0]) != -1 || index($data_annotation{$key}{'gene'}, $db1_cols[0]) != -1) {
                    $data_annotation{$key}{'other'} .= $db1_cols[0].";";
                }
            }
            unless (index($data_annotation{$key}{'transcript'}, $db1_cols[1]) != -1) {
                $data_annotation{$key}{'transcript'}.=$db1_cols[1].";";
            }
        }
        
        
        foreach my $db2_col (@db2_columns){
            if ($db2_col eq ".") {
                last;
            }
            my @db2_cols = split(/\|/,$db2_col);
            if (length($data_annotation{$key}{'gene'}) == 0) {
                $data_annotation{$key}{'gene'} = $db2_cols[0];
            }else {
                unless (index($data_annotation{$key}{'other'}, $db2_cols[0]) != -1 || index($data_annotation{$key}{'gene'}, $db2_cols[0]) != -1) {
                    $data_annotation{$key}{'other'} .= $db2_cols[0].";";
                }
            }
            
            unless (index($data_annotation{$key}{'transcript'}, $db2_cols[1]) != -1) {
                $data_annotation{$key}{'transcript'}.= $db2_cols[1].";";
            }
        }
        
        $data_annotation{$key}{'other'} =~s/;$//;
        $data_annotation{$key}{'transcript'}=~ s/;$//;
        
        if (length($data_annotation{$key}{'gene'}) < 1) {
            $data_annotation{$key}{'gene'} =".";
        }
        if (length($data_annotation{$key}{'transcript'}) < 1) {
            $data_annotation{$key}{'transcript'} =".";
        }
        if (length($data_annotation{$key}{'other'}) < 1) {
            $data_annotation{$key}{'other'} =".";
        }
        
    }
    close($fh);
    
    open(my $fho, ">$file.edit") or die $!;
    foreach my $key (keys %data_annotation){
        my @columns = split("_",$key);
        if (length($data_annotation{$key}{'rest'}) < 1) {
            print $fho "$columns[0]\t$columns[1]\t$columns[2]\t$data_annotation{$key}{'gene'}\t$data_annotation{$key}{'transcript'}\t$data_annotation{$key}{'other'}\t.\n";
        }else{
            print $fho "$columns[0]\t$columns[1]\t$columns[2]\t$data_annotation{$key}{'gene'}\t$data_annotation{$key}{'transcript'}\t$data_annotation{$key}{'other'}\t$data_annotation{$key}{'rest'}\n";
        }
    }
    close($fho);
    
    
    
    system("$script_dir/bin/bedtools sort -i $file.edit > $file");
    
    open($fhb,"<$file") or die $!;
    $line = <$fhb>;
    chomp($line);
    @tmp = split("\t",$line);
    $cols = scalar @tmp;
    close($fhb);
    
    my $header = "Chr Start Stop Gene RefSeq/HGMD Other_names";
    foreach my $i (7 .. $cols){
        my $tmp_cnt = $i-3;
        $header .= " Usr_col_$tmp_cnt";
    }
    
    system("echo $header > $index/HEADER.txt");
    
    system("awk \'{print > \"$index/\"\$1\".targetindex.txt\"}\' $file");
    my $rm_file = "$file.edit";
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = "$bed.Annotated-tmp";
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
}else{
    $file = "$default_index/vcrome2.1_hg19.flattened.bed";
    $index = $default_index;
}


foreach my $bam_file(@bam_files){
    my $tmp = $bam_file;
    #$tmp =~ s{.*/}{}; # remove path
    #$tmp = "$output_dir/$tmp";
    $tmp=~s/\.bam//;
    print STDERR "Generating cov.fasta file...\n\n";
    system("$covfasta_gen_script -o $tmp -t $file -i $bam_file");
    print STDERR "Running ExCID main script...\n\n";
    system("$script -f $tmp.cov.fasta -m $min -i $index");
}


### SUBROUTINES ###
sub USAGE {
    print "\nUSAGE: $0 -f <BAM file> -m <min threshold> -i <Bed file>\n";
    print "  -f:  BAM file. These option can be used multiple times for more than 1 bam\n";
    print "  -m:  minimal threshold; default 20\n";
    print "  -i:  Bed file\n\n";
    exit;
}