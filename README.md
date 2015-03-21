## About the ExCID Report ##

The Exome Coverage and Identification (ExCID) Report is a software tool developed at BCM-HGSC to assess sequence depth in user-defined targeted regions.  The tool was initially developed for use in targeted capture applications, but its functionality has evolved to encompass any sequencing application from amplicon and targeted capture sequencing to WGS.  ExCID analyzes sequence depth of any sequencing event, reports the average coverage across each target, and identifies bases below a user-defined threshold (20X coverage by default).  Furthermore, the tool annotates the target with the latest gene, transcript, and exon information from RefSeq and the Human Gene Mutation Database (HGMD).  The report has the option to output data tracks of sample targets and coverage that can be visualized in UCSC and IGV genome browsers.

## Outputs ##
* Outputs length, average coverage, and gene annotations for targets in BCM-HGSC VCRome (or your custom design)
* Outputs all regions of low coverage, including length and avg. coverage
* Output coverage track across regions of interest viewable in standard browser
* Outputs the percentage of gene covered in the design for all the Genetest genes (Clinically important genes) and other provided Gene Lists or Gene databases.

## Installation ##

Requirements:

        1. Latest version of JAVA and PERL.
        2. If on a Mac, you might need to install XCode: https://developer.apple.com/xcode/downloads/

1) Fill the information in Config.txt.

        DataBaseDir=/path/to/directory/to_put_the_databases/
        AnnotationDir=/path/to/directory/to_put_the_annotations_of_bed_files/
        
2) Run setup.sh script from command line.

        $./setup.sh
        
The setup script installs the bedtools version 2.17.0 (Released under GNU public license version 2 (GPL v2)) and maintained by the Quinlan Laboratory at the University of Virginia.
It will download the latest RefSEQ, VEGA, CCDS and miRBASE databses for bed file annotation. The databases generated are for Coding regions only.
        
## Usage ##

1) For using with VCrome regions run the program as:

        $ perl ExCID.BatchScript_v2.0-threading_Final.pl -bam <BAM file> -m <min threshold>
        
Multiple Bam files can be provided eg.
    
        $ perl ExCID.BatchScript_v2.0-threading_Final.pl -bam <BAM file1> -bam <BAM file2> -bam <BAM file3> -m <min threshold>
    OR
        $ perl ExCID.BatchScript_v2.0-threading_Final.pl -bamList <BAM list> -m <min threshold>
          where the BAM list is a test file with 1 bam file per line.
        
If the minimum threshold is not provided by the user then 20x coverage is assummed by default.
    
    
2) For using with a user defined BED file:

        $ perl ExCID.BatchScript_v2.0-threading_Final.pl -bam <BAM file> -m <min threshold> -i <Bed file>
        
Multiple Bam files can be provided eg.
    
        $ perl ExCID.BatchScript_v2.0-threading_Final.pl -bam <BAM file1> -bam <BAM file2> -bam <BAM file3> -m <min threshold> -i <Bed file>
        
If the minimum threshold is not provided by the user then 20x coverage is assummed by default.


3) For generating a wig file for all the target regions and a bed file for low covered regions for visualization in standard genome browser, use the '-wig' option:

        $ perl ExCID.BatchScript_v2.0-threading_Final.pl -bam <BAM file> -m <min threshold> -i <Bed file> -wig
        

4) Using '-d' option will not consider duplicates reads for generating the coverages statistics:

        $ perl ExCID.BatchScript_v2.0-threading_Final.pl -bam <BAM file> -m <min threshold> -i <Bed file> -d
        

The BED file will be annotated with RefSEQ, CCDS, VEGA and miRBASE gene annotations. The RefSEQ, CCDS, VEGA and miRBASE database can be updated as:
        $ ./update_databases.sh    


The Gentest genes were complied and annotated in November 2013.

## File Formats ##

1) If the user wants to generate a Gene database to obtain the Gene coverage percentage, the database should be of the following format:
    
    CHR START   STOP    GENE|TRASCRIPT_exon_number
    
    Example:
    10	100177320	100177483	HPS1|NM_000195_cds_0
    10	100177931	100178014	HPS1|NM_000195_cds_1
    10	100179801	100179915	HPS1|NM_000195_cds_2

