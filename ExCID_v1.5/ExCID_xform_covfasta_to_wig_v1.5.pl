#!/usr/bin/perl
#
# Description: create a wiggle file
#  from the cov.fasta file, so the
#  WGL can view the exome as a track
#  in UCSC.  Will create a bedfile 
#  from the lowcov annotated file for
#  use in UCSC.
#######################################
use strict;
use diagnostics;
use Getopt::Std;

### GLOBAL VARS ###
my %opt;
### OPTIONS ###
getopts("f:l:",\%opt);
my $infile = $opt{f} || &USAGE;
my $lowcov_in = $opt{l} || &USAGE; #meant for sites <20X coverage

#my $infile = $ARGV[0] || die "Please enter a file in cov.fasta format\n"; #input cap_stats.cov.fasta
#my $outfile = "cap_stats.VCRome.wig";
(my $outfile = $infile) =~ s{.*/}{}; # remove path
$outfile =~ s{\.[^.]+$}{}; # removes extension 
$outfile .= ".VCRome.wig";

# All output files will have this header:
open(OUTFILE, ">$outfile") || die "Can't open $outfile: $!\n";
print OUTFILE "track type=wiggle_0 name=Exome_Coverage autoScale=off viewLimits=0.0:100.0 graphType=bar visibility=full\n";


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

#my $lowcov_in = "cap_stats.lowcov-anno.txt";
#my $lowcov_out = "cap_stats.lowcov.bed";

(my $lowcov_out = $lowcov_in) =~ s{.*/}{}; # remove path
$lowcov_out =~ s{\.[^.]+$}{}; # removes extension 
$lowcov_out .= ".bed";

open(LOWOUT, ">$lowcov_out") || die "Can't open $lowcov_out: $!\n";
open(LOWIN, "<$lowcov_in") || die "Can't open $lowcov_in: $!\n";
print LOWOUT "track name=\'Low Coverage\' description=\'Coverage <20X\' color=255,0,0\n";
while (<LOWIN>) {
    chomp $_;
    if ($_ =~ m/^\#/) {next;}
    if ($_ =~ m/^[0-9]$/) {next;}
    my @low_array = split(/\t/, $_);
    my $lowchr = $low_array[0];
    my $lowstart = $low_array[1];
    my $lowstop = $low_array[2];
    print LOWOUT "$lowchr\t$lowstart\t$lowstop\n";
}
close LOWOUT;
close LOWIN;


sub USAGE {
    print "\nUSAGE: $0 -f <cov.fasta> -l <cap_stats.lowcov-anno.txt>\n";
    print "  -f:  file in fasta coverage format.\n";
    print "  -l:  cap_stats.lowcov-anno.txt generated from ExCID_getlowcov_and_annotate_exome.pl\n\n";
    exit;
}
