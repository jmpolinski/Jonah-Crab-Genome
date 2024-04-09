#!/bin/bash
##################################################################################################
##                        QC OF GENE MODEL PREDICTIONS TO GET FINAL SET                         ##
##         	       Code used by JMP -  Analysis performed in Sept 2023                      ##
##################################################################################################
cd /data/prj/JonahCrabGenome/Analyses/6_genePrediction/
mkdir qc && cd qc


# remove * stop codons from braker amino acid fasta
cp ../braker3_2023-07/braker.aa braker_noStop.faa
sed -i -e 's/[*]\+//g' braker_noStop.faa


# run interproscan to find predictions with known protein domains
/data/app/interproscan-5.63-95.0/interproscan.sh -appl Pfam -i braker_noStop.faa -iprlookup -goterms
awk '{print $1}' braker_noStop.faa.tsv | sort | uniq > keep_genes.w.pfam


# extract multi-exonic genes with GFACs
mkdir gfacs-multiex 
module load gFACs
gFACs.pl -f braker_2.0_gtf --fasta ../JC-genome_v3-230128.fasta.masked --rem-monoexonics --get-protein-fasta --get-fasta -O ./gfacs-multiex/ ../braker3_2023-07/braker.gtf

# extract multiexonic genes with no interproscan hit:
grep ">ID=g.*.\." genes.fasta.faa > multiex-genes_iso
grep ">" genes.fasta.faa > tmp
grep -vf multiex-genes_iso tmp > multiex-noIso
## need to add the t back in for transcript isoforms and .t1 for genes with no isoforms to make them match original fasta headers
sed -i -e 's/\./\.t/g' multiex-genes_iso
sed -i -e 's/$/.t1/g' multiex-noIso
cat multiex-noIso multiex-genes_iso | sort > multiex-genes
rm multiex-noIso multiex-genes_iso tmp
sed -i -e 's/>ID=//g' multiex-genes
cd ..
grep -vw -f ../keep_genes.w.pfam multiex-genes > multiex_noPfam
cd ../
module load seqtk
seqtk subseq ./braker_noStop.faa gfacs-multiex/multiex_noPfam > braker_multiexnoPfam.faa


# blast multiexonic genes without interproscan hit against SwissProt db and extract IDs for ones with hits:
module load blast
blastp -query braker_multiexnoPfam.faa -db /data/app/databases/uniprot/uniprot -num_threads 12 -outfmt 6 -evalue 1e-1 -out multi-noP_uniprot.txt
awk '{print $1}' multi-noP_uniprot.txt | sort | uniq > keep_multex-noP_SwissProt


# blast multiexonic genes without interproscan hit or SwissProt hit to NR database and extra IDs for ones with hits:
grep -wv -f keep_multex-noP_SwissProt gfacs-multiex/multiex_noPfam > multiex_noPfam_noSP
/data/app/seqtk/seqtk subseq braker_multiexnoPfam.faa multiex_noPfam_noSP > multiex_noPfam_noSP.faa
blastp -query ./multiex_noPfam_noSP.faa -db /data/app/blast/nr -num_threads 8 -outfmt 6 -evalue 1e-10 -out multi_noPnoSP_nr.txt -max_target_seqs 5
awk '{print $1}' multi_noPnoSP_nr.txt | sort | uniq > keep_mutliex-noP-w.NRhit


# Create a fasta file of gene subset with protein domain or that are multiexonic and have a blast hit (Swissprot + NR)
cat keep* > ../all_QC-genes.ID
/data/app/seqtk/seqtk subseq braker.aa all_QC-genes.ID > cbor_braker3-qc_isoforms.faa
/data/app/seqtk/seqtk subseq braker.codingseq all_QC-genes.ID > cbor_braker3-qc_isoforms.fa

# check BUSCO completeness
module load busco
busco -i cbor_braker3-qc_isoforms.faa -o BUSCO_qc.braker3 -l /data/app/busco-5.3.1/lineages/arthropoda_odb10 --mode protein -c 8

echo "Gene QC completed -- check busco results to see how it went!"
