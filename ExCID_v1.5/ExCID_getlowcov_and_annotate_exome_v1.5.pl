#!/usr/bin/perl -w
#
# Description: Given file with fasta.cov format,
#  determines avg coverage of each target region,
#  and detects bases lower than user-defined cov
#  threshhold.  Also determines avg-cov for the
#  poorly covered bases.
#
#  Annotation: additionally, for each target
#  and subsequent lowcoverage region, annotates
#  genename(s) and transcript-exons covered.
#
# Date: 2013 Feb 07; version1
# Date: 2013 April 20; version1.1 - NV: Added index as an input parameters, 
#   removed hardlinks, changed output filenames as input_file.VCRome-anno.txt 
# Date: 2013 May 24: version1.2 Bug Fixes
#   1) Fixed: targets not sitting on genes or mirna are skipped
#   2) Fixed: exact coverages of 20X being ignored.
# USAGE: $0 -f <fasta.cov> -m <min_cov> -i <path/to/capture_design_index_directory>
####################################################################################
use strict;
use diagnostics;
use Getopt::Std;

### GLOBAL VARS ###
my %opt;
my $chr;
my $start_coord;
my $stop_coord;
my $anno_result_out;
my $script_dir;
my $warn = 1;

### OPTIONS ###
getopts("f:m:i:",\%opt);
my $file = $opt{f} || &USAGE;
my $min = $opt{m} || 20; #meant for sites <20X coverage
my $index = $opt{i} || &USAGE;

### OUTPUT FILES ###
#my $VCRome_target_outfile = "cap_stats.VCRome-anno.txt";
my $VCRome_target_outfile = $file; 
#$VCRome_target_outfile =~ s{.*/}{}; # remove path
$VCRome_target_outfile =~ s{\.[^.]+$}{}; # removes extension 
$VCRome_target_outfile .= ".VCRome-anno.txt";

#my $VCRome_lowcov_outfile = "cap_stats.lowcov-anno.txt";
my $VCRome_lowcov_outfile = $file;
#$VCRome_lowcov_outfile =~ s{.*/}{}; # remove path
$VCRome_lowcov_outfile =~ s{\.[^.]+$}{}; # removes extension 
$VCRome_lowcov_outfile .= ".lowcov-anno.txt";

### BODY ###
my $script_dir_tmp = $0;
$script_dir_tmp =~m/^.+\//; 
$script_dir=$&;

open (HEADER, "$index/HEADER.txt") or die "Header file missing in $index:$!\n";
my $header = <HEADER>;
my @header_cols = split(" ",$header);
my $out_header = join("\t",@header_cols[0..2])."\tLength\tCoverage\t".join("\t",@header_cols[3..(scalar(@header_cols)-1)]);
open (OUTVCR, ">$VCRome_target_outfile") || die "Can't open $VCRome_target_outfile:$!\n";
open (OUTLOW, ">$VCRome_lowcov_outfile") || die "Can't open $VCRome_target_outfile:$!\n";
print OUTVCR '##generatedBy:', "$0", "\n";
print OUTVCR '##fileDate=', eval { `date +"%Y%m%d"` };   # not portable to non-Unix
print OUTVCR "##$out_header\n";
print OUTLOW '##generatedBy:', "$0", "\n";
print OUTLOW '##fileDate=', eval { `date +"%Y%m%d"` };   # not portable to non-Unix
print OUTLOW "##$out_header\n";

open (FILE, "<$file") || die "Can't open $file: $!\n";
while (<FILE>) {
    chomp $_;
    if ($_ =~ m/^\>/) {
        ($chr, $start_coord, $stop_coord) = split(/ /, $_);
        $chr =~ s/\>//;
	$chr =~ s/[Cc]hr//; #if the number has a 'chr' prefix, remove it.
        next;
    } elsif ($_ =~ m/^$/) {
	#this line is blank; skip it.
 	next; 
    }

    my @array = split(/ /, $_);
    my $array_size = scalar(@array);
    my $target_sum = eval join '+',@array;
    my $avg_target_cov = $target_sum / $array_size;
    $avg_target_cov = sprintf("%.2f", $avg_target_cov);
    #print OUTVCR "$chr\t$start_coord\t$stop_coord\t$array_size\t$avg_target_cov\n";
    $anno_result_out = &ANNO_TARGET($index, $chr, $start_coord, $stop_coord, $array_size, $avg_target_cov);
    print OUTVCR "$anno_result_out\n";
    my $array_marker = $start_coord + $array_size - 1;
    my $count = $start_coord;

    my $start;
    my $end;
    my $iterate;
    my $comment;

    my @lowcov_array;
    my $low_cov_array_size;
    my $low_cov_sum;
    my $low_cov_avg;

    foreach (@array) {
	    my $num = $_; chomp $_;
	    if ( $num >= $min ) 
	    {
	        if (( $iterate) && ( $iterate < $count )) {
	            $end = $count - 1;
		    $low_cov_array_size = @lowcov_array;
		    $low_cov_sum = eval join '+',@lowcov_array;
		    $low_cov_avg = $low_cov_sum / $low_cov_array_size;
		    $low_cov_avg = sprintf("%.2f", $low_cov_avg);
	            #print OUTLOW "$chr\t$start\t$end\t$low_cov_array_size\t$low_cov_avg\n";
		    $anno_result_out = &ANNO_TARGET($index, $chr, $start, $end, $low_cov_array_size, $low_cov_avg);
		    print OUTLOW "$anno_result_out\n";
	            $start = "";
	            $end = "";
	            $iterate = "";
	            $comment = "above threshold";
		    $low_cov_array_size = "";
		    $low_cov_sum = "";
		    $low_cov_avg = "";
	            $count++;
	        } else {
	            $count++;
	            $comment = "above threshold";
	        }
	    }
	    elsif ( $num < $min )
	    {

	        if ( ! $start ) {
	            $start = $count;
	            $iterate = $count;
		    push(@lowcov_array, $num);
	            $comment = "start & iteration started";
	        } elsif ( $count == $array_marker ) {
		    push(@lowcov_array, $num);
		    $low_cov_array_size = @lowcov_array;
		    $low_cov_sum = eval join '+',@lowcov_array; 
		    $low_cov_avg = $low_cov_sum / $low_cov_array_size;
		    $low_cov_avg = sprintf("%.2f", $low_cov_avg);
	            #print OUTLOW "$chr\t$start\t$count\t$low_cov_array_size\t$low_cov_avg\n";
		    $anno_result_out = &ANNO_TARGET($index, $chr, $start, $count, $low_cov_array_size, $low_cov_avg);
		    print OUTLOW "$anno_result_out\n";
	        } elsif ( $start <  $count ) {
	            $iterate++;
	            $comment = "iterating";
		    push(@lowcov_array, $num);
	        }
	        $count++;
	    }
	    #print "$num\t$count\t$comment\n";
	}
}
close FILE;
close OUTLOW;
close HEADER;

### SUBROUTINES ###
sub USAGE {
    print "\nUSAGE: $0 -f <cov.fasta> -m <min threshold> -i </path/to/capture_design_index_directory>\n";
    print "  -f:  file in fasta coverage format.\n";
    print "  -m:  minimal threshold; default 20\n";
    print "  -i:  path to capture design index directory\n\n";
    exit;
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
	        $annoresult = "$lowchr\t$lowstart\t$lowstop\t$length\t$cov";
		for(my $i = 3; $i<scalar(@columns);$i++){
		    $annoresult .="\t$columns[$i]";
		}
		return $annoresult;
	    }
    }
    if (! $annoresult) {
	$annoresult = "$lowchr\t$lowstart\t$lowstop\t$length\t$cov";
	return $annoresult;
    }
    close INDEX;
}
