#!/bin/bash

BASEDIR=$(cd `dirname ${0}`; pwd)
source Config.txt

mkdir -p $DataBaseDir; 
mkdir -p $AnnotationDir;
rm -rf $DataBaseDir/* ;
rm -rf $AnnotationDir/* ;


rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz $DataBaseDir ;
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ccdsGene.txt.gz $DataBaseDir ;
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/vegaGene.txt.gz $DataBaseDir ;
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/vegaGtp.txt.gz $DataBaseDir ;
cp miRBASE_r20.gff2 $DataBaseDir/miRBASE_r20.gff2 ;

ls --color=never $DataBaseDir/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
ls --color=never $DataBaseDir/*.txt | while read FILE ; do awk -F "\t" '{print $3"\t"$7"\t"$8"\t"$2"\t"$13"\t"$9"\t"$10"\t"$11}' "$FILE" > "$FILE.bed" ; done  ;
rm $DataBaseDir/vegaGtp.txt.bed ;
awk -F "\t" '{print $1"\t"$2}' $DataBaseDir/vegaGtp.txt > $DataBaseDir/VEGA-hgnc_names ; 
rm $DataBaseDir/vegaGtp.txt ;

awk -F "\t| " '{print $1"\t"$4"\t"$5"\t"$10}' $DataBaseDir/miRBASE_r20.gff2  | sed s/ID=\"//   | sed s/\"\;//  | grep "^#" -v > $DataBaseDir/miRBASE_r20.bed ;

perl creat_bed_UCSC_coding.pl $DataBaseDir/refGene.txt.bed ;
perl creat_bed_UCSC_coding.pl $DataBaseDir/ccdsGene.txt.bed ;
perl creat_bed_UCSC_coding.pl $DataBaseDir/vegaGene.txt.bed ;

awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$5}' $DataBaseDir/ccdsGene.txt.bed-exon.bed > tmp ; mv tmp $DataBaseDir/ccdsGene.txt.bed-exon.bed ;

sed s/^chr// $DataBaseDir/miRBASE_r20.bed > tmp ; mv tmp $DataBaseDir/miRBASE_r20.bed ;

perl Get_HGNC.pl  > $DataBaseDir/HGNC_database.txt ;

perl check_HGNC_individual_mirnaDB.pl $DataBaseDir/miRBASE_r20.bed $DataBaseDir/HGNC_database.txt > $DataBaseDir/miRBASE_r20_HGNC.bed  &
perl check_HGNC_individual_CCDSDB.pl $DataBaseDir/ccdsGene.txt.bed-exon.bed $DataBaseDir/HGNC_database.txt > $DataBaseDir/ccdsGene.txt.bed-exon_HGNC.bed  &
perl check_HGNC_individual_VEGADB.pl $DataBaseDir/vegaGene.txt.bed-exon.bed $DataBaseDir/HGNC_database.txt $DataBaseDir/VEGA-hgnc_names > $DataBaseDir/vegaGene.txt.bed-exon_HGNC.bed  &
perl check_HGNC_individual_RefSeqDB.pl $DataBaseDir/refGene.txt.bed-exon.bed $DataBaseDir/HGNC_database.txt > $DataBaseDir/refGene.txt.bed-exon_HGNC.bed &

wait;

#grep -P "\tNM_" $DataBaseDir/refGene.txt.bed-exon_HGNC.bed | $BASEDIR/bin/bedtools intersect -a - -b $DataBaseDir/refGene.txt.bed-Coding_region.bed -u > refGene.txt.bed-exon_HGNC.bed_tmp ;
#grep -P "\tNM_" $DataBaseDir/refGene.txt.bed-exon_HGNC.bed -v | cat - refGene.txt.bed-exon_HGNC.bed_tmp > tmp ;
#mv tmp $DataBaseDir/refGene.txt.bed-exon_HGNC.bed ;

#$BASEDIR/bin/bedtools intersect -a $DataBaseDir/vegaGene.txt.bed-exon_HGNC.bed -b $DataBaseDir/vegaGene.txt.bed-Coding_region.bed -u > tmp;
#mv tmp $DataBaseDir/vegaGene.txt.bed-exon_HGNC.bed ;


#rm refGene.txt.bed-exon_HGNC.bed_tmp ;
rm $DataBaseDir/miRBASE_r20.gff2 ;
rm $DataBaseDir/miRBASE_r20.bed ;
rm $DataBaseDir/ccdsGene.txt ;
rm $DataBaseDir/vegaGene.txt ;
rm $DataBaseDir/refGene.txt ;
rm $DataBaseDir/ccdsGene.txt.bed ;
rm $DataBaseDir/vegaGene.txt.bed ;
rm $DataBaseDir/refGene.txt.bed ;
rm $DataBaseDir/ccdsGene.txt.bed-exon.bed ;
rm $DataBaseDir/vegaGene.txt.bed-exon.bed ;
rm $DataBaseDir/refGene.txt.bed-exon.bed ;

