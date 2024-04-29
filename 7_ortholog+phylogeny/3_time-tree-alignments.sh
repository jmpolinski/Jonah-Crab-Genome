#!/bin/bash
################## ALIGNMENTS FOR TIME-SCALED PHYLOGENETIC TREE ##################

cd /data/prj/Fisheries/JonahCrabGenome/Analyses/10_orthology/timetree

# fastas with single-copy orthologs from OrthoFinder
ln -s ../OrthoFinder/Results_Nov08_2/Single_Copy_Orthologue_Sequences/ .
ls Single_Copy_Orthologue_Sequences > files

# alignemnts
module load muscle/v5.1
for i in `cat files`; do muscle -align Single_Copy_Orthologue_Sequences/$i -output muscle-aligned-SCO/${i}.afa -threads 16; done

sed -i -e 's/.fa//g' files

# making sure that all alignments have species in the same order
mkdir aligned-sorted-SCO
module load seqkit/v2.5.1
for in in `cat files`; do seqkit sort -i ./muscle-aligned-SCO/${i}.fa.afa > aligned-sorted-SCO/${i}_sorted.afa; done

# create list of alignment files for Gblock input
ls aligned-sorted-SCO/ > aligned
sed -i -s 's$OG$aligned-sorted-SCO/OG$g' aligned

# run Gblocks
module load gblocks/v0.91b
Gblocks aligned -t=p -v=10000 -a=y
mv aligned-sorted-SCO/aligned-gb.seq ./aligned_all-gb.afa

echo "Export aligned_all-gb.afa and use in MEGA"
