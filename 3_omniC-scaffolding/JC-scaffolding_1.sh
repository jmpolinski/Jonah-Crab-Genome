#!/bin/bash
##################################################################################################
##                         Scaffolding with Omni-C data (part 1)				##
##	   Code used by JMP to generate a chrosome-level assembly for C. borealis		##
## 			     Analysis performed in August 2022					##
##################################################################################################

cd /data/prj/Fisheries/JonahCrabGenome/Analyses/2_OmniC/juicer

## Step 1. Remove contigs <10 Kb from preliminary assembly
perl scripts/removesmall.pl 10000 JC_flye-v2.fasta > JC-draft_10kb+.fa

## Step 2. Juicer analysis
# Generate BWA index for draft genome:
bwa index JC-draft_10kb+.fa

# Generate restriction site file:
/data/app/juicer-1.6/misc/generate_site_positions.py none JC-draft_10kb+ JC-draft_10kb+.fa

# Generate chromosome size file:
samtools faidx JC-draft_10kb+.fa
cut -f 1,2 JC-draft_10kb+.fa.fai > chrom.sizes

#Run Juicer:
/data/app/juicer-1.6/CPU/juicer.sh -g JC-draft_10kb+ -s "none" -z JC-draft_10kb+.fa -p chrom.sizes -t 24


## Step 3. 3D-DNA analysis
/data/app/3d-dna-201008/run-asm-pipeline.sh JC-draft_10kb+.fa ./aligned/merged_nodups.txt

## Step 4. Visualize and modify in JBAT
echo "The next step is done in JBAT manually. Transfer and open the files JC-draft_10kb+.final.hic and JC-draft_10kb+_HiC.assembly in JBAT." 
echo "Check the chromosome designation in JBAT, manually correct as needed, and export the modified assembly as JC-draft_10kb+.final.review.assembly."
echo "Once that is complete, proceed to the next script -- JC-scaffolding_2.sh"

