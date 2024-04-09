# Gene Model Prediction

Prior to ab initio gene model prediction with BRAKER, RNA-seq data from 4 tissues and 3 individuals were mapped to soft-masked genome.  
RNAseq data was location in the directory /data/prj/JonahCrabGenome/Data/RNAseq/.  

## RNA-seq mapping with Hisat2
First built a Hisat2 index for the reference genome:
```
ln -s ../JC-genome_v3-230128.fasta.masked
hisat2-build -f JC-genome_v3-230128.fasta JC-genome_v3
```

A text files names "samples" was created that had the IDs for the three individuals used for RNAseq:
```
JC502
JC504
JC506
```
  
A bash script ```1_rnaseq-hisat2.sh``` was written t unzip fastqs, map reads to the genome, gzip fastqs, and convert SAM files to BAM files. It was written so each tissue type could be mapped simultaneously in loops:
```
# gill 
for i in `cat samples`; do ./1_rnaseq-hisat2.sh $i Gill /data/prj/JonahCrabGenome/Data/RNAseq/; done
# hemocytes
for i in `cat samples`; do ./1_rnaseq-hisat2.sh $i Hem /data/prj/JonahCrabGenome/Data/RNAseq/; done
# hepatopancreas
for i in `cat samples`; do ./1_rnaseq-hisat2.sh $i Hepato /data/prj/JonahCrabGenome/Data/RNAseq/; done
# muscle
for i in `cat samples`; do ./1_rnaseq-hisat2.sh $i Muscle /data/prj/JonahCrabGenome/Data/RNAseq/; done
```
Individual BAMs were combined into a single BAM to use a evidence in BRAKER:
```
module load samtools
samtools merge -@ 32 JC_all-RNAseq.bam JC50*.bam
```

## BRAKER gene model prediction
Initial gene model prediction was done with BRAKER2 (v2.1.6) in January 2023. 
However, it was noted in July 2023 that a new version of BRAKER had been released (BRAKER3 - v3.0.3). 
Ab initio gene model prediction was rerun with BRAKER3 using the hintsfile.gff from the initial run with BRAKER2. 
Both runs used Augustus v3.4.0. 
  
  
### BRAKER2 preliminary run
```
module load braker/v2.1.6
braker.pl --genome=JC-genome_v3-230128.fasta.masked --bam=JC_all-RNAseq.bam --softmasking --cores=16 --species=Cborealis_braker-rna --workingdir=rna-hints
module unload braker/v2.1.6
```
### BRAKER3 re-run:
```
module load braker/v3.0.3
braker.pl --genome=JC-genome_v3-230128.fasta.masked --hints=./rna-hints/hintsfile.gff --threads=16 --species=Cborealis_braker-rna --useexisting 
```
  
  
## Gene model QC

A script was written to first run interproscan and save gene model IDs for those with known protein domains. 
Additional genes were kept based on the criteria of being multiexonic and having high homology (e<1e-10) to proteins in either the SwissProt or NR databases.
```
./4_gene-model-qc.sh
```

  
## Functional Annotation and Isoform selection
The newest version of RefSeq was downloaded 9/8/2023 and used for EnTAP annotation
```
# first make diamond database for refseq
cd /data/resources/databases/refseq_2023-07
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete.nonredundant_protein.{1..445}.protein.faa.gz
cat complete*protein.faa.gz > refseq-completeNR-protein_2023.faa.gz
rm complete*
gunzip refseq-completeNR-protein_2023.faa.gz
module load diamond/v2.1.8
diamond makedb --in refseq-completeNR-protein_2023.faa -d refseq-protein_completeNR-2023
gzip refseq-completeNR-protein_2023.faa
ln -s /data/resources/databases/refseq_2023-07/refseq-protein_completeNR-2023.dmnd /data/app/EnTAP-0.10.8-beta/databases/refseq-prot.dmnd  

# run entap
module load entap/v0.10.8
cd /data/prj/Fisheries/JonahCrabGenome/Analyses/7_geneAnnotation
EnTAP --runP -i ../6_genePrediction/qc/cbor_braker3-qc_isoforms.faa -d /data/resources/databases/refseq_2023-07/refseq-protein_completeNR-2023.dmnd \
 --ini /data/resources/app_modules/EnTAP-v0.10.8-beta/entap_config.ini -t 16

```
For genes with multiple isoforms, the isoform with the highest homology to a RefSeq protein was selected. If no isoforms for a gene had a RefSeq match, the longest isoform was kept.
  

## Final results: 
The final gene model set included 24,830 gene models and exhibited a BUSCO score of 92.3% [Single-copy = 88.1%, Duplicated = 4.2%]
