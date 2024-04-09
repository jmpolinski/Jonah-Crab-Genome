#!/bin/bash
##################################################################################################
##                         MAPPING RNASEQ DATA TO THE REFERENCE                                 ##
##         Code used by JMP to generate evidence for ab initio gene prediction                  ##
##                           Analysis performed in January 2023                                 ##
##################################################################################################


# set a tissue type
SAM=$1
TIS=$2
LOC=$3
if [[ -z ${SAM} && -z ${TIS} && -z ${LOC} ]];
then echo "Missing input to run command. Please provide an TISSUE TYPE to put at front of output files"; exit;
else echo "Variables provided:" ${SAM} ${TIS} ${LOC}; fi

# check current directory
pwd

# set import variables/path
export HISAT2_INDEXES=/data/prj/JonahCrabGenome/RNAseq/hisat2-mapping/
echo "HISAT2_INDEXES" $HISAT2_INDEXES
export PATH=/data/app/hisat2-2.2.1/:$PATH
echo "PATH" $PATH
pwd

# unzip reads 
echo "unzipping fastqs:" ${SAM} ${TIS}
date
gunzip ${LOC}/${SAM}-${TIS}*.gz

# begin mapping loop
echo "mapping with hisat2:" ${SAM} ${TIS}
date
/data/app/hisat2-2.2.1/hisat2 -p 8 -x JC-genome_v3 -1 ${LOC}/${SAM}-${TIS}*_R1_001.fastq -2 ${LOC}/${SAM}-${TIS}*_R2_001.fastq -S ${SAM}-${TIS}.sam

# compress fastq files after mapping
echo "gzipping fastqs:" ${SAM} ${TIS}
date
gzip ${LOC}/${SAM}-${TIS}*.fastq

# convert sam to bam & removing sam
echo "converting SAM to BAM:" ${SAM} ${TIS}
date
samtools view -Su ${SAM}-${TIS}.sam | samtools sort -o ${SAM}-${TIS}.bam
rm ${SAM}-${TIS}.sam

#complete!
echo "mapping and cleanup complete for:" ${SAM} ${TIS}
date

