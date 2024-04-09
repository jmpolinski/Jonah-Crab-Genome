#!/bin/bash
##################################################################################################
##                     PRELIMINARY AB INITIO PREDICTION WITH BRAKER2                            ##
##         Code used by JMP to generate evidence for ab initio gene prediction                  ##
##                           Analysis performed in January 2023                                 ##
##################################################################################################

module load braker/v2.1.6
braker.pl --genome=JC-genome_v3-230128.fasta.masked --bam=JC_all-RNAseq.bam --softmasking --cores=16 --species=Cborealis_braker-rna --workingdir=rna-hints

