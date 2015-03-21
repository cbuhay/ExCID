#!/usr/bin/perl -w
## Anootater for BED file.
use strict;
use diagnostics;
use Getopt::Std;
use File::Basename;

### GLOBAL VARS ###
my $bed = $ARGV[0] || &USAGE;
my $database = $ARGV[1] || &USAGE;
my $output_dir = $ARGV[2] || dirname($bed);
my $db_name = basename($database);

system("ln -s $bed $output_dir/");

$bed = "$output_dir/".basename($bed);
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
$cmd .= "\$".($cols+3+1)."\"\t\""."\$".($cols+3+1+1)."\"\t\""."\$".($cols+3+1+1+1)." ";


system("bedtools intersect -a $bed -b $database -wao | awk -F\$\'\t\' '{print $cmd}' > $bed.$db_name.Annotated");

my %data_annotation;
open(my $fh, "< $bed.$db_name.Annotated") or die $!;
while (my $line = <$fh>) {
    chomp($line);
    $line=~s/\t-1\t/\t.\t/;
    my @columns = split("\t",$line);
    my $key = "$columns[0]_$columns[1]_$columns[2]";
    unless (exists $data_annotation{$key}){
        $data_annotation{$key}{'HGNC_gene'} = "";
        $data_annotation{$key}{'transcript'} = "";
        $data_annotation{$key}{'other_Genes'} = "";
        $data_annotation{$key}{'rest'} = join("\t", @columns[3..(scalar(@columns)-4)]);
    }
    
    if ($columns[scalar(@columns)-1] ne ".") {
        if (length($data_annotation{$key}{'HGNC_gene'}) == 0) {
            $data_annotation{$key}{'HGNC_gene'} = $columns[scalar(@columns)-1].";";
        }else{
            my $check = $columns[scalar(@columns)-1].";";
            if (index($data_annotation{$key}{'HGNC_gene'},$check) == -1 && index($data_annotation{$key}{'other_Genes'},$check) == -1) {
                $data_annotation{$key}{'other_Genes'} = $columns[scalar(@columns)-1].";";
            }
        }
    }
    
    if ($columns[scalar(@columns)-2] ne ".") {
        if (length($data_annotation{$key}{'transcript'}) == 0) {
            $data_annotation{$key}{'transcript'} = $columns[scalar(@columns)-2].";";
        }else{
            my $check = $columns[scalar(@columns)-2].";";
            unless (index($data_annotation{$key}{'transcript'}, $check) != -1) {
                $data_annotation{$key}{'transcript'} .= $columns[scalar(@columns)-2].";";
            }
        }
    }
    
    
    if ($columns[scalar(@columns)-3] ne ".") {        
        my $check = $columns[scalar(@columns)-3].";";
        if (index($data_annotation{$key}{'HGNC_gene'},$check) == -1 && index($data_annotation{$key}{'other_Genes'},$check) == -1) {
            $data_annotation{$key}{'other_Genes'} = $columns[scalar(@columns)-3].";";
        }
    }
    
    
    if (length($data_annotation{$key}{'HGNC_gene'}) == 0) {
        $data_annotation{$key}{'HGNC_gene'} =".";
    }
    if (length($data_annotation{$key}{'transcript'}) == 0) {
        $data_annotation{$key}{'transcript'} =".";
    }
    if (length($data_annotation{$key}{'other_Genes'}) == 0) {
        $data_annotation{$key}{'other_Genes'} =".";
    }      
}
close($fh);

open(my $fho, ">$bed.$db_name.Annotated.edit") or die $!;
foreach my $key (keys %data_annotation){
    my @columns = split("_",$key);
    
    $data_annotation{$key}{'HGNC_gene'} =~ s/;$//;
    $data_annotation{$key}{'other_Genes'} =~ s/;$//;
    $data_annotation{$key}{'transcript'}=~ s/;$//;
    $data_annotation{$key}{'other_Genes'} =~ s/^\.// if ($data_annotation{$key}{'other_Genes'} ne ".");
    
    if (length($data_annotation{$key}{'rest'}) < 1) {
        my $out_line = "$columns[0]\t$columns[1]\t$columns[2]\t$data_annotation{$key}{'HGNC_gene'}\t$data_annotation{$key}{'transcript'}\t$data_annotation{$key}{'other_Genes'}\t.\n";
        $out_line=~s/\t-1\t/\t.\t/;
        print $fho $out_line;
    }else{
        my $out_line = "$columns[0]\t$columns[1]\t$columns[2]\t$data_annotation{$key}{'HGNC_gene'}\t$data_annotation{$key}{'transcript'}\t$data_annotation{$key}{'other_Genes'}\t$data_annotation{$key}{'rest'}\n";
        $out_line=~s/\t-1\t/\t.\t/;
        print $fho $out_line;
    }
}
close($fho);



system("bedtools sort -i $bed.$db_name.Annotated.edit > $bed-$db_name");

my $rm_file = "$bed.$db_name.Annotated.edit";
unlink $rm_file or warn "Could not unlink $rm_file: $!";
$rm_file = "$bed.$db_name.Annotated";
unlink $rm_file or warn "Could not unlink $rm_file: $!";


### SUBROUTINES ###
sub USAGE {
    print "USAGE: $0 <Bed file> <Database file> <output dir>\n\n";
    exit;
}