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
my $database1 = "$script_dir/database/refGene.txt.bed-exon_HGNC.bed";
my $database2 = "$script_dir/database/ccdsGene.txt.bed-exon_HGNC.bed";
my $database3 = "$script_dir/database/vegaGene.txt.bed-exon_HGNC.bed";
my $database4 = "$script_dir/database/miRBASE_r20_HGNC.bed";
my $combined='';


### OPTIONS ###
GetOptions('i:s' => \$opt{i}, 'l:s' => \$opt{l}, 'list:s' => \$opt{g}, 'db:s' => \$opt{db}) || &USAGE;
my $bed = $opt{i} || "null";
my $low_cov_file = $opt{l} || "null";
my $gene_list = $opt{g} || "null";
my $targeted_database = $opt{db} || "null" ;
my $pct_file;
my @genes_lest;


if ($gene_list eq "null" && $targeted_database eq "null") {
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




if($targeted_database && -e $targeted_database ){
    print STDERR  "Targeted Data base provided. Obtaining Gene percent coverages.\n";
    $pct_file = genes_check_db($low_cov_file, $bed, $targeted_database);
    exit 0;
}else{
    print STDERR "Some issue with the run. Please check the commands.\n";
    exit -1;
}



########## SUBROUTINES ###########
sub USAGE {
    print "\nUSAGE: $0 -i <Targeted BED file> -l <low cov file> -list <Genelist> -db <Gene database>\n";
    print "  -i:  Annotated bed file\n";
    print "  -l:  Low cov bed file\n";
    print "  -list:  List of genes. one per line\n";
    print "  -db:  Database of interested Genes. Please see documentation for database format. \n";
    exit;
}

sub genes_check_db {
    
    my ($low_cov_file, $anno_bed, $db_list)  = @_;
    my $db_name = basename($db_list);
    my $out_file = $low_cov_file;
    $out_file=~ s{.*/}{}; # remove path
    $out_file.="-miRBASE-anno.bed";
    my $nottargeted_rfs = $anno_bed;
    $nottargeted_rfs=~ s{.*/}{}; # remove path
    $nottargeted_rfs.="-notTargeted_in_DB-$out_file.bed";
    my $tmp = $out_file;
    $tmp=~ s{.*/}{}; # remove path
    $tmp=~s/\.bed$//;
    $combined = $nottargeted_rfs;
    $combined=~s/\-$out_file//;
    $combined=~s/\.bed$//;
    $combined.="-$tmp.bed";
    
    system("$script_dir/bin/bedtools subtract -a $db_list -b $anno_bed > $nottargeted_rfs");
    system("$script_dir/bin/bedtools intersect -a $db_list -b $low_cov_file > $out_file");
    system("cat $nottargeted_rfs $out_file | $script_dir/bin/bedtools sort -i  > $combined");
    system("$script_dir/bin/bedtools merge -i $combined > $combined-tmp");
    system("$script_dir/bin/bedtools intersect -a $db_list -b $combined-tmp > $combined");
    my $db_list_transcript_size = get_Transcirpt_size_db($db_list);
    my $out_FILE = get_Transcirpt_size($low_cov_file,$combined, $db_list_transcript_size);
    my $rm_file = "$combined-tmp";
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = $out_file;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = $nottargeted_rfs;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    $rm_file = "$combined";
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    
    $db_name=~ s{\.[^.]+$}{}; # removes extension
    system("mv $out_FILE $db_name-pct.txt");
    return "$db_name-pct.txt";
    
    return $out_FILE;
}
    
sub get_Transcirpt_size{
    my ($low_cov_file,$combined,$db)=@_;
    
    my $tmp = $combined;
    $tmp=~ s{.*/}{}; # remove path
    $tmp=~s/\.bed$//;
    
    
    my $outfile_final = "$tmp"."_transcriptSIZE.txt";
    open(my $in,"$combined") || die "Can't open $combined: $!\n";
    open(my $out, ">$outfile_final") || die "Can't open $outfile_final: $!\n";
    
    while(<$in>)
         {
             chomp; my $line = $_; 
             my ($chr, $start, $stop, $gene) = split(/\s/, $line);
             my $exonsize = $stop - $start + 1;
             my @tranarray = split(/\;/, $gene);
             my $arraysize = scalar(@tranarray) - 1;
             for (my $j=0; $j<=$arraysize; $j++) {
                print $out "$tranarray[$j]\t$exonsize\t$chr\t$start\t$stop\n";
             }
         }
    
    close($in);
    close($out);
    my $outFile = get_pct($outfile_final,$db);
    return $outFile;
}

sub get_Transcirpt_size_db{
    my ($combined)=@_;
    
    my $tmp = $combined;
    $tmp=~ s{.*/}{}; # remove path
    $tmp=~s/\.bed$//;
    
    
    my $outfile_final = "$tmp"."_transcriptSIZE.txt";
    open(my $in,"$combined") || die "Can't open $combined: $!\n";
    open(my $out, ">$outfile_final") || die "Can't open $outfile_final: $!\n";
    
    while(<$in>)
         {
             chomp; my $line = $_; 
             my ($chr, $start, $stop, $gene) = split(/\s/, $line);
             my $exonsize = $stop - $start + 1;
             my @tranarray = split(/\;/, $gene);
             my $arraysize = scalar(@tranarray) - 1;
             for (my $j=0; $j<=$arraysize; $j++) {
                print $out "$tranarray[$j]\t$exonsize\t$chr\t$start\t$stop\n";
             }
         }
    
    close($in);
    close($out);
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
        if (scalar @data1 ==1) {
            $control{"$data1[0]_$data1[2]_$data1[3]"}{"val"} = 0;
            $control{"$data1[0]_$data1[2]_$data1[3]"}{"Gene"} = $data[0];
            $control{"$data1[0]_$data1[2]_$data1[3]"}{"ID"} = $data[1];
        }
        if (scalar @data1 ==5) {
            $control{"$data1[0]_$data1[2]_$data1[3]"}{"val"} = $data1[1];
            $control{"$data1[0]_$data1[2]_$data1[3]"}{"Gene"} = $data[0];
            $control{"$data1[0]_$data1[2]_$data1[3]"}{"ID"} = $data[1];
            $control{"$data1[0]_$data1[2]_$data1[3]"}{"chr"} = $data1[2];
            $control{"$data1[0]_$data1[2]_$data1[3]"}{"start"} = $data1[3];
            $control{"$data1[0]_$data1[2]_$data1[3]"}{"stop"} = $data1[4];
        }
    }
    close($fh);
    
    open(my $fh1,"<$infile") or die $!;
    while (my $line = <$fh1>) {
        chomp($line);
        my @data1 = split("\t",$line);
        my @data = split(/\|/,$data1[0]);
        if (scalar @data1 ==1) {
            $file_data{"$data1[0]_$data1[2]_$data1[3]"}{"val"} = 0;
            $file_data{"$data1[0]_$data1[2]_$data1[3]"}{"Gene"} = $data[0];
            $file_data{"$data1[0]_$data1[2]_$data1[3]"}{"ID"} = $data[1];
        }
        if (scalar @data1 ==5) {
            $file_data{"$data1[0]_$data1[2]_$data1[3]"}{"val"} = $data1[1];
            $file_data{"$data1[0]_$data1[2]_$data1[3]"}{"Gene"} = $data[0];
            $file_data{"$data1[0]_$data1[2]_$data1[3]"}{"ID"} = $data[1];
        }
    }
    close($fh1);
    
    my @low_region_keys = keys %file_data;
    open(my $fho,">$outfile") or die $!;
    
    foreach my $key (keys %control){
        my @key_split = split(/_/,$key);
        my $no_match = 1;
        
        my @matches = grep { /$key_split[0]/ } @low_region_keys;
        foreach my $match (@matches){
            my @low_key_split = split(/_/,$match);
            
            if (($key_split[0] eq $low_key_split[0]) && ($key_split[1] eq $low_key_split[1]) && $key_split[2] <= $low_key_split[2]) {
                my $tmp = $file_data{$match}{"val"}/$control{$key}{"val"};
                my $pct = sprintf("%.3f",(1-$tmp)*100);
                
                my $grep_key = $control{$key}{"Gene"}."|".$control{$key}{"ID"};
                my @regions = `grep -w \"$grep_key\" $combined `;
                my $lowcov_coords='';
                foreach my $lowcov_region (@regions){
                    my @tmp = split ("\t",$lowcov_region);
                    if ("$tmp[0]" eq $control{$key}{"chr"} && $tmp[1]>= $control{$key}{"start"} && $tmp[2] <= $control{$key}{"stop"}) {
                        $lowcov_coords.="$tmp[1]-$tmp[2];";
                    }
                }
                $lowcov_coords=~s/;$//;
                
                print $fho $control{$key}{"chr"}."\t".$control{$key}{"Gene"}."\t".$control{$key}{"ID"}."\t1\t".$control{$key}{"start"}."\t".$control{$key}{"stop"}."\t$pct%\t$lowcov_coords\n";
                $no_match = 0;
            }
            
        }
        if ($no_match == 1) {
            print $fho $control{$key}{"chr"}."\t".$control{$key}{"Gene"}."\t".$control{$key}{"ID"}."\t1\t".$control{$key}{"start"}."\t".$control{$key}{"stop"}."\t100.000%\t.\n";
        }
    }
    close($fho);
    my $rm_file = $infile;
    unlink $rm_file or warn "Could not unlink $rm_file: $!";
    
    return $outfile;
}

##############################################################################################################################

sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d", $t->year + 1900, $t->mon + 1, $t->mday, $t->hour, $t->min, $t->sec );
}