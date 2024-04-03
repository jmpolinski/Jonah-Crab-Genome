# Removing mitochondrial genome contamination from the nuclear genome
  
## Step 1. Find mtDNA in genome
In order to remove any contamination originating from mitochrondrial genome sequence that was misassembled into the nuclear genome, two references from the same genus were used:
- [*Cancer magister*](https://doi.org/10.1080/23802359.2019.1691474) - GenBank:[MN371144.1](https://www.ncbi.nlm.nih.gov/nuccore/MN371144.1)
- [*Cancer pagurus*](https://doi.org/10.1080/23802359.2019.1689859) - GenBank: [MN334534.1](https://www.ncbi.nlm.nih.gov/nuccore/MN334534.1)  
  
Both reference mtDNA sequences were combined into a file ```references.fa```. This was used for a blast search (1_blast-mt.sh)
  
## Step 2. Assess where mtDNA was assembled
Unfortunately, the mtDNA hits appear to be misassembled within larger scaffolds. As the reference sequences has already been shared with colleagues performing lc-WGS, it was determined that the best way forward was to hard-mask (N) the mitochondiral contamination. 

## Step 3. Hard-masking mtDNA
Working in: ```JonahCrabGenome/Analyses/3_mitoGenome/2_corrected-scafs```

1) Link needed files:
```
ln -s ../../JC-genome_v1_mitoG.fasta
ln -s ../../JC-genome_v1_mitoG.fasta.fai
```
2) Created a text files with the chromosome and coordinates to be replaced with N's (used multiples of 80 since the genome files lines are 80 characters)
3) for loop to pull out those regions:
```
for i in `cat mito-regions-to-N`; do samtools faidx JC-genome_v1_mitoG.fasta ${i} > ${i}_N.fa; done
for i in `cat regions-to-keep`; do samtools faidx JC-genome_v1_mitoG.fasta ${i} > ${i}_ok.fa; done
```
4) substitute A/C/T/G's for N's:
```
sed -i -e 's/A/N/g' *_N.fa
sed -i -e 's/T/N/g' *_N.fa
sed -i -e 's/C/N/g' *_N.fa
sed -i -e 's/G/N/g' *_N.fa
```
5) recombine contigs and then remove extra headers:
```
cat Cbor_chr20:* > Cbor_chr20-corrected.fa
nano Cbor_chr20-corrected.fa
```
6) switch " | " in original genome file to "--" so HiC scaffold info isn't lost, then separate out chromosomes into separate fastas:
```
sed -i -e 's/ | HiC/--HiC/g' JC-genome_v1_mitoG.fasta
grep ">" JC-genome_v1_mitoG.fasta > scaffolds
sed -i -e 's/>//g' scaffolds
nano scaffolds ## remove scaffolds 2, 20, 31, 36, 41
samtools faidx JC-genome_v1_mitoG.fasta
for i in `cat scaffolds`; do samtools faidx JC-genome_v1_mitoG.fasta $i > $i.fa; done
```
7) Recombine all scaffolds and move to main directory:
```
cat *fa > JC-genome_v3-230128.fasta
mv JC-genome_v3-230128.fasta* ../../
```
