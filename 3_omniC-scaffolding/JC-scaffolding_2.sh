#!/bin/bash
##################################################################################################
##                         Scaffolding with Omni-C data (part 2)				##
##	   Code used by JMP to generate a chrosome-level assembly for C. borealis		##
## 			     Analysis performed in August 2022					##
##################################################################################################

cd /data/prj/Fisheries/JonahCrabGenome/Analyses/2_OmniC/juicer/3d-dna/juicebox-reviewed

## Step 5. Output final revised assembly with 3D-DNA using reviewed assembly file from JBAT
/data/app/3d-dna-201008/run-asm-pipeline-post-review.sh -r ./JC-draft_10kb+.final.review.assembly ../../JC-draft_10kb+.fa ../../aligned/merged_nodups.txt
