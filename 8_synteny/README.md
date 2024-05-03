# Synteny analysis 

Three species were used for synteny analyses: Jonah crab (*Cancer borealis*), Chinese mitten crab (*Eriocheir sinensis*), and swimming crab (*Portunus trituberculatus*). 
These were chosen for inclusion based on being the crab species with BUSCO scores >85% complete and <10% duplicated. 
  
[MCscan (python version)](https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)) was used to conduct this analysis. 

  
## Step 1 - prepare input files
pwd = ```/data/prj/Fisheries/JonahCrabGenome/Analyses/8_synteny```
### Jonah crab:
```
module load jcvi_MCscan/v1.3.8
python -m jcvi.formats.gff bed --type=transcript ../cbor_braker_v2.gff -o cborealis.bed
python -m jcvi.formats.fasta format ../cbor_braker3-qc_v2.fa cborealis.cds

#subset out only genes on the 51 chromosome-length scaffolds
grep "Cbor_chr" cborealis.bed > cborealis_chr.bed
awk '{print $4}' cborealis_chr.bed > cbor_chr-only-genes
module load seqtk
seqtk subseq cborealis.cds cbor_chr-only-genes > cborealis_chr.fa
python -m jcvi.formats.fasta format cborealis_chr.fa cborealis_chr.cds
```
### Chinese mitten crab:
```
mkdir Esinensis && cd Esinensis/
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_024679095.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_024679095.1.zip" -H "Accept: application/zip"
unzip *zip
mv ./*/*/* .
module load jcvi_MCscan/v1.3.8
python -m jcvi.formats.gff bed --type=mRNA --key=Name Esinensis/genomic.gff -o esinensis.bed
python -m jcvi.formats.fasta format Esinensis/rna.fna esinensis.cds
# remove all but ID in CDS file headers
sed -i -e 's/\s.*$//g' esinensis.cds

# chromosome scaffolds only
grep ">*chromosome " Esinensis/GCF_024679095.1_ASM2467909v1_genomic.fna | sed 's/>//g' | sed 's/ Eriocheir*.*//g' > Esinensis_chrIDs
grep -f Esinensis_chrIDs esinensis.bed > esinensis_chr.bed
awk '{print $4}' esinensis_chr.bed > esin_chr-only-genes
module load seqtk
seqtk subseq esinensis.cds esin_chr-only-genes > esinensis_chr.fa
python -m jcvi.formats.fasta format esinensis_chr.fa esinensis_chr.cds
```
### Swimming crab:
```
mkdir ptrituberculatus && cd ptrituberculatus
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_017591435.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_017591435.1.zip" -H "Accept: application/zip"
unzip GCF_017591435.1.zip
mv  ncbi_dataset/data/GCF_017591435.1/* .
cd ..

# gff to bed for MCscan
module load jcvi_MCscan/v1.3.8
python -m jcvi.formats.gff bed --type=mRNA --key=Name ptrituberculatus/genomic.gff -o ptrituberculatus.bed

# get chromosome-level scaffolds only
grep ">*chromosome " ptrituberculatus/GCF_017591435.1_ASM1759143v1_genomic.fna | sed 's/>//g' | sed 's/ Portunus*.*//g' > Ptrituberculatus_chrIDs
grep -f Ptrituberculatus_chrIDs ptrituberculatus.bed > ptrituberculatus_chr.bed
awk '{print $4}' ptrituberculatus_chr.bed > ptri_chr-only-genes
module load seqtk
seqtk subseq ptrituberculatus/rna.fna ptri_chr-only-genes > ptrituberculatus_chr.fa
sed -i -e 's/ PRED.*.//g' ptrituberculatus_chr.fa
python -m jcvi.formats.fasta format ptrituberculatus_chr.fa ptrituberculatus_chr.cds

python -m jcvi.formats.fasta format Esinensis/rna.fna esinensis.cds
# remove all but ID in CDS file headers
sed -i -e 's/\s.*$//g' esinensis.cds
```

  
  
## Pairwise synteny: Jonah crab + Chinese mitten crab 
```
# creates orthology  files and oxford plots
python -m jcvi.compara.catalog ortholog cborealis_chr esinensis_chr --no_strip_names

# synteny depth histograms
python -m jcvi.compara.synteny depth --histogram cborealis_chr.esinensis_chr.anchors

# synteny karyotypes
python -m jcvi.compara.synteny screen --minspan=10 --simple cborealis_chr.esinensis_chr.anchors cborealis_chr.esinensis_chr.anchors.new
python -m jcvi.graphics.karyotype seqids_c+e layout_c+e --chrstyle=roundrect
mv karyotype.pdf cborealis_chr.esinensis_chr.karyotype.pdf
```

## Pairwise synteny: Jonah crab + Swimming crab 
```
# creates orthology files and oxford plots
python -m jcvi.compara.catalog ortholog ptrituberculatus_chr cborealis_chr  --no_strip_names

# synteny depth histograms
python -m jcvi.compara.synteny depth --histogram ptrituberculatus_chr.cborealis_chr.anchors

# synteny karyotypes
python -m jcvi.compara.synteny screen --minspan=10 --simple ptrituberculatus_chr.cborealis_chr.anchors ptrituberculatus_chr.cborealis_chr.anchors.new
python -m jcvi.graphics.karyotype seqids_p+c layout_p+c --chrstyle=roundrect
mv karyotype.pdf ptrituberculatus_chr.cborealis_chr.karyotype.pdf
```

## Pairwise synteny: Swimming crab + Chinese mitten crab
```
# creates orthology files and oxford plots
python -m jcvi.compara.catalog ortholog ptrituberculatus_chr esinensis_chr  --no_strip_names

# synteny depth histograms
python -m jcvi.compara.synteny depth --histogram ptrituberculatus_esinensis_chr.anchors

# synteny karyotypes
python -m jcvi.compara.synteny screen --minspan=10 --simple ptrituberculatus_chr.esinensis_chr.anchors ptrituberculatus_chr.esinensis_chr.anchors.new
python -m jcvi.graphics.karyotype seqids_p+e layout_p+e --chrstyle=roundrect
mv karyotype.pdf ptrituberculatus_chr.esinensis_chr.karyotype.pdf
```
