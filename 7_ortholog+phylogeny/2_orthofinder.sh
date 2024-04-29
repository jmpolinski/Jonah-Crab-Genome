#!/bin/bash
################### ORTHOFINDER ANALYSIS ###################
#        run on 11/7/2023 using Orthofinder v2.5.4         #
############################################################

cd /data/prj/Fisheries/JonahCrabGenome/Analyses/10_orthology

# Fastas with protein sequences for species included in this analysis are location in /data/prj/Fisheries/JonahCrabGenome/Analyses/8_gainLoss/proteins

module load orthofinder/v2.5.4
orthofinder -t 16 -a 16  -f ./proteins
