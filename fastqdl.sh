#!/bin/bash

study_folder=$1 # The order was fixed in the main wgs workflow script.
accession=$2

echo "DOWNLOADING DATA FROM : ${accession} , INTO : ${study_folder}/00.rawdata/"

echo "fastq-dl START"
source activate fastq-dl

# the ${study_folder} variable is calling the entire pathway under the hood
fastq-dl --accession ${accession} \
	 --outdir "${study_folder}/00.rawdata" \
	 --cpus 20 \

conda deactivate
echo "fastq-dl FINISHED"
