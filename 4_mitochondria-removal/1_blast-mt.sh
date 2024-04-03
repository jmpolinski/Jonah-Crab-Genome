#!/bin/bash

#################################################
#           blast search for mtDNA              #
#################################################

cd /data/prj/Fisheries/JonahCrabGenome/Analyes/3_mitoGenome/

module load blast
blastn -query references.fa -subject ../JC-genome_all.fasta -outfmt 6 -out jc_mitoGenome-blastn.txt -num_threads 8 -evalue 1e-10
