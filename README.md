## About the ExCID Report ##

The Exome Coverage and Identification (ExCID) Report has been developed to analyze the sequence depth of each exome, report the average coverage across each target, and identify bases below 20X.  Furthermore, the tool also reports the gene, transcript, and exon captured by the target.  The report has the option to output data tracks of sample targets and coverage that can be visualized in UCSC and IGV genome browsers.

## Outputs ##
* Outputs length, average coverage, and gene annotations for targets in BCM-HGSC VCRome (or your custom design)
* Outputs all regions of low coverage, including length and avg. coverage
* Output coverage track across regions of interest viewable in standard browser

## Installation ##

Requirements: Latest version of JAVA and PERL.

1) Run setup.sh script from command line.
        $./setup.sh
        
The setup script installs the bedtools version 2.17.0 (Released under GNU public license version 2 (GPL v2)) and maintained by the Quinlan Laboratory at the University of Virginia. 
        
## Usage ##

1) For using with VCrome regions run the program as:

        $ perl ExCID.BatchScript.pl -f <BAM file> -m <min threshold>
        
Multiple Bam files can be provided eg.
    
        $ perl ExCID.BatchScript.pl -f <BAM file1> -f <BAM file2> -f <BAM file3> -m <min threshold>
        
If the minimum threshold is not provided by the user then 20x coverage is assummed by default.
    
    
2) For using with a user defined BED file:

        $ perl ExCID.BatchScript.pl -f <BAM file> -m <min threshold> -i <Bed file>
        
Multiple Bam files can be provided eg.
    
        $ perl ExCID.BatchScript.pl -f <BAM file1> -f <BAM file2> -f <BAM file3> -m <min threshold> -i <Bed file>
        
If the minimum threshold is not provided by the user then 20x coverage is assummed by default.
The BED file will be annotated with RefSEQ and HGMD gene annotations. The RefSEQ and HGMD database was obtained in November 2013.

To update the RefSEQ and HGMD database please look at the files in database directory for the database format. These are all the entries in the respective databases.
