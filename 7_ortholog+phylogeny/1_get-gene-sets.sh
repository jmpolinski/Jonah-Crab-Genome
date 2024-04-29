#!/bin/bash 
# Selection of Species for Ortholog Analysis

cd /data/prj/Fisheries/JonahCrabGenome/Analyses/10_orthology

mkdir proteins && cd proteins
module load busco/v5.4.5
module load gFACs/v1.1.2


### Jonah crab 
ln -s ../../Cbor_geneModels_v2.faa ./Cborealis_proteins.faa
busco -i Cborealis_proteins.faa -o ./lobster_busco -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:92.3%[S:88.1%,D:4.2%],F:1.8%,M:5.9%,n:1013
sed -i -e 's/>/>Cborealis_/g' Cborealis_proteins.faa


### [American lobster (*Homarus americanus*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018991925.1/)
mkdir hamericanus && cd hamericanus
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_018991925.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_018991925.1.zip" -H "Accept: application/zip"
unzip GCF_018991925.1.zip
sed -i -e 's/ Homarus.*.//g' ncbi_dataset/data/GCF_018991925.1/GCF_018991925.1_GMGI_Hamer_2.0_genomic.fna
gFACs.pl -p Hamericanus --fasta ncbi_dataset/data/GCF_018991925.1/GCF_018991925.1_GMGI_Hamer_2.0_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_018991925.1/genomic.gff
busco -i genes.fasta.faa -o hamericanus-refseq -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:96.4%[S:95.4%,D:1.0%],F:2.1%,M:1.5%,n:1013 
sed -i -e 's/ID=gene-/Hamericanus_/g' genes.fasta.faa
mv genes.fasta.faa ../Hamericanus_proteins.faa


### [Swimming crab (*Portunus trituberculatus*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_017591435.1/)
mkdir ptrituberculatus && cd ptrituberculatus
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_017591435.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_017591435.1.zip" -H "Accept: application/zip"
unzip GCF_017591435.1.zip
sed -i -e 's/ Portunus.*.//g' ncbi_dataset/data/GCF_017591435.1/GCF_017591435.1_ASM1759143v1_genomic.fna
gFACs.pl -p ptrituberculatus--fasta ncbi_dataset/data/GCF_017591435.1/GCF_017591435.1_ASM1759143v1_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_017591435.1/genomic.gff
busco -i genes.fasta.faa -o ptrituberculatus -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:97.0%[S:96.4%,D:0.6%],F:1.4%,M:1.6%,n:1013 
sed -i -e 's/ID=gene-/Ptrituberculatus_/g' genes.fasta.faa
mv genes.fasta.faa ../Ptrituberculatus_proteins.faa


### [Pacific white shrimp (*Litopenaeus vannamei*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003789085.1/)
mkdir lvannamei && cd lvannamei
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_003789085.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_003789085.1.zip" -H "Accept: application/zip"
unzip GCF_003789085.1.zip
sed -i -e 's/ Penaeus.*.//g' ncbi_dataset/data/GCF_003789085.1/GCF_003789085.1_ASM378908v1_genomic.fna
gFACs.pl -p Lvannamei --fasta ncbi_dataset/data/GCF_003789085.1/GCF_003789085.1_ASM378908v1_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_003789085.1/genomic.gff
busco -i genes.fasta.faa -o lvannamei-refseq -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:86.2%[S:77.9%,D:8.3%],F:5.3%,M:8.5%,n:1013
sed -i -e 's/ID=gene-/Lvannamei_/g' genes.fasta.faa
mv genes.fasta.faa ../Lvannamei_proteins.faa


### [Mud crab (*Scylla paramamosain*)](https://figshare.com/articles/dataset/Supportting_files_for_the_article_of_A_chromosome-level_genome_of_the_mud_crab_Scylla_paramamosain_Estampador_provides_insights_into_the_evolution_of_chemical_and_light_perception_in_this_crustacean_published_in_Molecular_Ecology_Resources_/13338968)
# downloaded file "SP-protein.fa" from FigShare
busco -i SP-protein.fa -o SP-protein -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:72.0%[S:66.6%,D:5.4%],F:2.7%,M:25.3%,n:1013


### [Water flea (*Daphnia magna*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_020631705.1/)  
mkdir dmagna && cd dmagna 
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_020631705.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_020631705.1.zip" -H "Accept: application/zip"
unzip GCF_020631705.1.zip
# need to remove all but scaffold ID from head for gFACs
sed -i -e 's/ Daphnia.*.//g' ncbi_dataset/data/GCF_020631705.1/GCF_020631705.1_ASM2063170v1.1_genomic.fna
gFACs.pl -p Dmagna --fasta ncbi_dataset/data/GCF_020631705.1/GCF_020631705.1_ASM2063170v1.1_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_020631705.1/genomic.gff
busco -i genes.fasta.faa -o dmagna-gFACS_busco -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10/ --mode protein -c 8
   # C:97.9%[S:92.0%,D:5.9%],F:0.2%,M:1.9%,n:1013
sed -i -e 's/>ID=/>Dmagna_/g' Dmagna_proteins.faa
mv genes.fasta.faa ../Dmagna_proteins.faa


### [Fruit fly (*Drosophila melanogaster*)](https://ftp.flybase.net/releases/current/dmel_r6.54/gff/)  
mkdir dmelanogaster && cd dmelanogaster 
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001215.4/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000001215.4.zip" -H "Accept: application/zip"
unzip GCF_000001215.4.zip
# need to remove all but scaffold ID for gFACs to work
sed -i -e 's/ Dro.*.//g' ncbi_dataset/data/GCF_000001215.4/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
gFACs.pl -p Dmelanogaster --fasta ncbi_dataset/data/GCF_000001215.4/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna -f refseq_gff --unique-genes-only --statistics --get-protein-fasta proteins.faa ./ -O ./ ncbi_dataset/data/GCF_000001215.4/genomic.gff
busco -i genes.fasta.faa -o busco-gfacs -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:99.7%[S:98.9%,D:0.8%],F:0.1%,M:0.2%,n:1013
sed -i -e 's/ID=gene-//g' Dmelanogaster_proteins.faa
mv genes.fasta.faa ../Dmelanogaster_proteins.faa


### [Amphipod crustacean (*Hyalella azteca)*](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000764305.2/)
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000764305.2/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000764305.2.zip" -H "Accept: application/zip"
unzip GCF_000764305.2.zip
# need to remove all but scaffold ID for gFACs to work
sed -i -e 's/ Hyal.*.//g' ncbi_dataset/data/GCF_000764305.2/GCF_000764305.2_Hazt_2.0.2_genomic.fna
gFACs.pl -p Hazteca --fasta ncbi_dataset/data/GCF_000764305.2/GCF_000764305.2_Hazt_2.0.2_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_000764305.2/genomic.gff
busco -i genes.fasta.faa -o hazteca-gfacs_busco -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10/ --mode protein -c 8
   # C:92.9%[S:92.2%,D:0.7%],F:2.6%,M:4.5%,n:1013
sed -i -e 's/ID=gene-/Hazteca_/g' Hazteca_proteins.faa
mv genes.fasta.faa ../Hazteca_proteins.faa


### [Water flea (*Daphnia pulex*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_021134715.1/)
mkdir dpulex && cd dpulex
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_021134715.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_021134715.1.zip" -H "Accept: application/zip"
unzip GCF_021134715.1.zip
head ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna
gFACs.pl -p Dpulex --fasta ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_021134715.1/genomic.gff
busco -i genes.fasta.faa -o dpulex-gfacs -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10/ --mode protein -c 8
   # C:97.4%[S:96.5%,D:0.9%],F:0.6%,M:2.0%,n:1013
sed -i -e 's/ID=gene-/Dpulex_/g' genes.fasta.faa
mv genes.fasta.faa ../Dpulex_proteins.faa


### [Chinese mitten crab (*Eriocheir sinensis*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_024679095.1/)
mkdir Esinensis && cd Esinensis
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_024679095.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_024679095.1.zip" -H "Accept: application/zip"
unzip GCF_024679095.1.zip
sed -i -e 's/ Eriocheir.*.//g' ncbi_dataset/data/GCF_024679095.1/GCF_024679095.1_ASM2467909v1_genomic.fna
gFACs.pl -p Esinensis--fasta ncbi_dataset/data/GCF_024679095.1/GCF_024679095.1_ASM2467909v1_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_024679095.1/genomic.gff
busco -i genes.fasta.faa -o esinensis-gfacs -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10/ --mode protein -c 8
   # C:97.5%[S:91.2%,D:6.3%],F:0.4%,M:2.1%,n:1013
sed -i -e 's/ID=gene-/Esinensis_/g' genes.fasta.faa
mv genes.fasta.faa ../Esinensis_proteins.faa
mv esinensis-gfacs/short_summary.specific.arthropoda_odb10.esinensis-gfacs.txt ../busco-summaries/Esinensis_short_summary.specific.arthropoda_odb10.txt


### [Black tiger shrimp (*Penaeus monodon*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_015228065.2/)
mkdir pmonodon && cd pmonodon
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_015228065.2/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_015228065.2.zip" -H "Accept: application/zip"
unzip GCF_015228065.2.zip
sed -i -e 's/ Penaeus.*.//g' ncbi_dataset/data/GCF_015228065.2/GCF_015228065.2_NSTDA_Pmon_1_genomic.fna
gFACs.pl -p Pmonodon --fasta ncbi_dataset/data/GCF_015228065.2/GCF_015228065.2_NSTDA_Pmon_1_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_015228065.2/genomic.gff
busco -i genes.fasta.faa -o pmonodon -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:80.1%[S:77.5%,D:2.6%],F:5.5%,M:14.4%,n:1013
sed -i -e 's/ID=gene-/Pmonodon_/g' genes.fasta.faa
mv genes.fasta.faa ../Pmonodon_proteins.faa


### [Red swamp crayfish (*Procambarus clarkii*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_020424385.1/)
mkdir pclarkii && cd pclarkii
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_020424385.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_020424385.1.zip" -H "Accept: application/zip"
unzip GCF_020424385.1.zip
sed -i -e 's/ Procambarus.*.//g' ncbi_dataset/data/GCF_020424385.1/GCF_020424385.1_ASM2042438v2_genomic.fna
gFACs.pl -p Pclarkii --fasta ncbi_dataset/data/GCF_020424385.1/GCF_020424385.1_ASM2042438v2_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_020424385.1/genomic.gff
busco -i genes.fasta.faa -o pclarkii -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10/ --mode protein -c 8
   # C:97.3%[S:95.3%,D:2.0%],F:1.4%,M:1.3%,n:1013
sed -i -e 's/ID=gene-/Pclarkii_/g' genes.fasta.faa
mv genes.fasta.faa ../Pclarkii_proteins.faa


### [Australian red claw crayfish (*Cherax quadricarinatus*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_026875155.1/)
mkdir cquadricarinatus && cd cquadricarinatus
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_026875155.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_026875155.1.zip" -H "Accept: application/zip"
unzip GCF_026875155.1.zip
sed -i -e 's/ Cherax.*.//g' ncbi_dataset/data/GCF_026875155.1/GCF_026875155.1_ASM2687515v2_genomic.fna
gFACs.pl -p Cquadricarinatus --fasta ncbi_dataset/data/GCF_026875155.1/GCF_026875155.1_ASM2687515v2_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_026875155.1/genomic.gff
busco -i genes.fasta.faa -o cquadricarinatus -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:83.7%[S:82.3%,D:1.4%],F:7.4%,M:8.9%,n:1013
sed -i -e 's/ ID=gene-/Cquadricarinatus_/g' genes.fasta.faa
mv genes.fasta.faa ../Cquadricarinatus_proteins.faa


### [Tidepool copepod (*Tigriopus californicus*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_007210705.1/)
mkdir tcalifornicus && cd tcalifornicus
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_007210705.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_007210705.1.zip" -H "Accept: application/zip"
unzip GCF_007210705.1.zip
sed -i -e 's/ Tigriopus.*.//g' ncbi_dataset/data/GCF_007210705.1/GCF_007210705.1_Tcal_SD_v2.1_genomic.fna
gFACs.pl -p tcalifornicus --fasta ncbi_dataset/data/GCF_007210705.1/GCF_007210705.1_Tcal_SD_v2.1_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_007210705.1/genomic.gff
busco -i genes.fasta.faa -o tcalifornicus -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10/ --mode protein -c 8
   # C:94.2%[S:92.7%,D:1.5%],F:1.1%,M:4.7%,n:1013
sed -i -e 's/ID=gene-/Tcalifornicus_/g' genes.fasta.faa
mv genes.fasta.faa ../Tcalifornicus_proteins.faa


### [Salmon louse (*Lepeophtheirus salmonis*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016086655.3/)
mkdir lsalmonis && cd lsalmonis
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_016086655.3/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_016086655.3.zip" -H "Accept: application/zip"
unzip GCF_016086655.3.zip
sed -i -e 's/ Lepeophtheirus.*.//g' ncbi_dataset/data/GCF_016086655.3/GCF_016086655.3_UVic_Lsal_1.2_genomic.fna
gFACs.pl -p lsalmonis --fasta ncbi_dataset/data/GCF_016086655.3/GCF_016086655.3_UVic_Lsal_1.2_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_016086655.3/genomic.gff
busco -i genes.fasta.faa -o lsalmonis -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:96.0%[S:92.5%,D:3.5%],F:0.9%,M:3.1%,n:1013 
sed -i -e 's/ID=gene-/Lsalmonis_/g' genes.fasta.faa
mv genes.fasta.faa ../Lsalmonis_proteins.faa


### [Brine shrimp (*Artemia franciscana*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032884065.1/*)
mkdir afranciscana && cd afranciscana
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_032884065.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCA_032884065.1.zip" -H "Accept: application/zip"
unzip GCA_032884065.1.zip
sed -i -e 's/ Artemia.*.//g' ncbi_dataset/data/GCA_032884065.1/GCA_032884065.1_ASM3288406v1_genomic.fna
gFACs.pl -p afranciscana --fasta ncbi_dataset/data/GCA_032884065.1/GCA_032884065.1_ASM3288406v1_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCA_032884065.1/genomic.gff
busco -i genes.fasta.faa -o afranciscana -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:81.0%[S:75.1%,D:5.9%],F:4.9%,M:14.1%,n:1013
sed -i -e 's/ID=gene-/Afranciscana_/g' genes.fasta.faa
mv genes.fasta.faa ../Afranciscana_proteins.faa


### [Fleshy prawn (*Penaeus chinensis*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_019202785.1/)
mkdir pchinensis && cd pchinensis
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_019202785.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_019202785.1.zip" -H "Accept: application/zip"
unzip GCF_019202785.1.zip
sed -i -s 's/ Penaeus.*.//g' ncbi_dataset/data/GCF_019202785.1/GCF_019202785.1_ASM1920278v2_genomic.fna
gFACs.pl -p pchinensis --fasta ncbi_dataset/data/GCF_019202785.1/GCF_019202785.1_ASM1920278v2_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_019202785.1/genomic.gff
busco -i genes.fasta.faa -o pchinensis -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:94.9%[S:93.8%,D:1.1%],F:1.4%,M:3.7%,n:1013
sed -i -e 's/ID=gene-/Pchinensis_/g' genes.fasta.faa
mv genes.fasta.faa ../Pchinensis_proteins.faa


### [Kuruma shrimp (*Penaeus japonicus*)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_017312705.1/)
mkdir pjaponicus && cd pjaponicus
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_017312705.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_017312705.1.zip" -H "Accept: application/zip"
unzip GCF_017312705.1.zip
sed -i -e 's/ Penaeus.*.//g' ncbi_dataset/data/GCF_017312705.1/GCF_017312705.1_Mj_TUMSAT_v1.0_genomic.fna
sed -i -e 's/ Marsupenaeus.*.//g' ncbi_dataset/data/GCF_017312705.1/GCF_017312705.1_Mj_TUMSAT_v1.0_genomic.fna
gFACs.pl -p pjaponicus --fasta ncbi_dataset/data/GCF_017312705.1/GCF_017312705.1_Mj_TUMSAT_v1.0_genomic.fna -f refseq_gff --unique-genes-only --get-protein-fasta test.faa ./ -O ./ ncbi_dataset/data/GCF_017312705.1/genomic.gff
busco -i genes.fasta.faa -o pjaponicus -l /data/resources/databases/BUSCO-lineages/arthropoda_odb10 --mode protein -c 8
   # C:95.5%[S:92.4%,D:3.1%],F:1.7%,M:2.8%,n:1013
sed -i -e 's/ID=gene-/Pjaponicus_/g' genes.fasta.faa
mv genes.fasta.faa ../Pjaponicus_proteins.faa

