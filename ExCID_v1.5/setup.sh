#!/bin/bash

BASEDIR=$(cd `dirname ${0}`; pwd)
echo "Do you wish to install BEDTools program?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) tar xzvf external_programs/BEDTools.v2.17.0.tar.gz 1>>setup.log 2>>setup.log ; 
              # build the bedtools package
              cd bedtools-2.17.0 ;
              make 1>>setup.log 2>>setup.log ;
              cp bin/bedtools ../bin/ ;
              break;;
        No )  echo "Please enter the path to the installed BEDTools bin directory (/Library/bedtools-2.17.0/bin/): " ;
              read input_variable ;
              cp $input_variable/bedtools $BASEDIR/bin/ ;
              break;;
    esac
done

cd $BASEDIR;
printf "\nExtracting VCRome indices... ";
tar -xzf $BASEDIR/VCRome_index.tgz;
printf "DONE\n";
printf "Extracting external databases... ";
tar -xzf $BASEDIR/database.tgz;
printf "DONE.\n\nSetup complete.\n";
