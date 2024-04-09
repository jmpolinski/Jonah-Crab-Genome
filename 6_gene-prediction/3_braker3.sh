#!/bin/bash
##################################################################################################
##                        RE-RUN AB INITIO PREDICTION WITH BRAKER3                              ##
##         Code used by JMP to generate evidence for ab initio gene prediction                  ##
##                           Analysis performed in July 2023	                                ##
##################################################################################################

module load braker/v3.0.3
braker.pl --genome=JC-genome_v3-230128.fasta.masked --hints=./rna-hints/hintsfile.gff --threads=16 --species=Cborealis_braker-rna --useexisting --workingdir=braker3_2023-07

