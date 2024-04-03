# Scaffolding with Omni-C data
  
Omni-C data was generated by Dovetail Genomics (Cantata Bio, LLC) and was used to scaffold the preliminary assembly generated by Flye.  
  
Juicer (v1.6), 3D-DNA (201008), and JuiceBox Assembly Tools (JBAT) were used to perform this analysis.  
  
## Step 1. Remove contigs <10 Kb from preliminary assembly
```
perl scripts/removesmall.pl 10000 JC_flye-v2.fasta > JC-draft_10kb+.fa
```
  
## Step 2. Juicer analysis  
Generate BWA index for draft genome:
```
bwa index JC-draft_10kb+.fa
```
  
Generate restriction site file:
```
/data/app/juicer-1.6/misc/generate_site_positions.py none JC-draft_10kb+ JC-draft_10kb+.fa
```
  
Generate chromosome size file:
```
samtools faidx JC-draft_10kb+.fa
cut -f 1,2 JC-draft_10kb+.fa.fai > chrom.sizes
```
  
Run Juicer:
```
/data/app/juicer-1.6/CPU/juicer.sh -g JC-draft_10kb+ -s "none" -z JC-draft_10kb+.fa -p chrom.sizes -t 24
```
  
## Step 3. 3D-DNA analysis
```
/data/app/3d-dna-201008/run-asm-pipeline.sh JC-draft_10kb+.fa merged_nodups.txt
```
  
## Step 4. Visualize and modify in JBAT
The HiC map showed ~51 distinct boxes, but 3D-DNA grouped all as a single chromosome. Manual review involved splitting the assembly up into the 51 segments with JBAT, exporting the new assembly, and getting new fasta, .assembly, and .hic files based on the corrected versions.  
  
## Step 5. Output final revised assembly with 3D-DNA using reviewed assembly file from JBAT
```
/data/app/3d-dna-201008/run-asm-pipeline-post-review.sh -r ./JC-draft_10kb+.final.review.assembly ./JC-draft_10kb+.fa ./merged_nodups.txt
```
  
  
  
The above described code can also be found in the two scripts "JC-scaffolding_1.sh" and "JC-scaffolding_2.sh". To use these scripts, users will need to adjust program/file locations. 