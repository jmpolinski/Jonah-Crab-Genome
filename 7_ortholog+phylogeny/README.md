# Ortholog and Phylogenetic Analyses

## Choosing species to include in analyses
[NCBI GenBank](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=6657) was searched for Crustacean genomes with NCBI RefSeq annotations. 
RefSeq annotations were downloaded as genome and GFF files. gFACs was used to extract gene models with only the longest isoform. BUSCO was used to check gene set completeness. 
Code for this can be found in ```1_get-gene-sets.sh```.  
  
BUSCO results (included species):  
| Name                | Species                    | Busco results                                | Classification				|
|:-------------------:|:--------------------------:|:--------------------------------------------:|:-------------------------------------------:|
| Jonah crab          | *Cancer borealis*          | C:92.3%[S:88.1%,D:4.2%],F:1.8%,M:5.9%,n:1013 | Arthropoda:Crustacea:Malacostraca:Decapoda  |
|American Lobster     | *Homarus americanus*       | C:96.4%[S:95.4%,D:1.0%],F:2.1%,M:1.5%,n:1013 | Arthropoda:Crustacea:Malacostraca:Decapoda  |
|Swimming crab        | *Portunus trituberculatus* | C:97.0%[S:96.4%,D:0.6%],F:1.4%,M:1.6%,n:1013 | Arthropoda:Crustacea:Malacostraca:Decapoda  |
|Pacific White shrimp | *Litopenaeus vannamei*     | C:86.2%[S:77.9%,D:8.3%],F:5.3%,M:8.5%,n:1013 | Arthropoda:Crustacea:Malacostraca:Decapoda  |
|Chinese Mitten crab  | *Eriocheir sinensis*       | C:97.5%[S:91.2%,D:6.3%],F:0.4%,M:2.1%,n:1013 | Arthropoda:Crustacea:Malacostraca:Decapoda  |
|Red swamp crayfish   | *Procambarus clarkii*      | C:97.3%[S:95.3%,D:2.0%],F:1.4%,M:1.3%,n:1013 | Arthropoda:Crustacea:Malacostraca:Decapoda  |
|Fleshy prawn         | *Penaeus chinensis*        | C:94.9%[S:93.8%,D:1.1%],F:1.4%,M:3.7%,n:1013 | Arthropoda:Crustacea:Malacostraca:Decapoda  |
|Kuruma shrimp        | *Penaeus japonicus*        | C:95.5%[S:92.4%,D:3.1%],F:1.7%,M:2.8%,n:1013 | Arthropoda:Crustacea:Malacostraca:Decapoda  |
|Salmon louse         | *Lepeophtheirus salmonis*  | C:96.0%[S:92.5%,D:3.5%],F:0.9%,M:3.1%,n:1013 | Arthropoda:Crustacea:Copepoda               |
|Tidepool copepod     | *Tigriopus californicus*   | C:94.2%[S:92.7%,D:1.5%],F:1.1%,M:4.7%,n:1013 | Arthropoda:Crustacea:Copepoda               |
|amphipod             | *Hyalella azteca*          | C:92.9%[S:92.2%,D:0.7%],F:2.6%,M:4.5%,n:1013 | Arthropoda:Crustacea:Malacostraca:Amphipoda |
|water flea           | *Daphnia magna*            | C:97.9%[S:92.0%,D:5.9%],F:0.2%,M:1.9%,n:1013 | Arthropoda:Crustacea:Branchiopoda           |
|water flea           | *Daphnia pulex*            | C:97.4%[S:96.5%,D:0.9%],F:0.6%,M:2.0%,n:1013 | Arthropoda:Crustacea:Branchiopoda           |
|Fruit fly            | *Drosophila melanogaster*  | C:99.7%[S:98.9%,D:0.8%],F:0.1%,M:0.2%,n:1013 | Arthropoda:Hexapoda:Insecta (outgroup)      |
|Bumblebee            | *Bombus pyrosoma*          | C:99.5%[S:99.5%,D:0.0%],F:0.0%,M:0.5%,n:1013 | Arthropoda:Hexapoda:Insecta (outgroup)      |




BUSCO results (species not included; <85% complete):
| Name              | Species                  | Busco results                                 |
|:-----------------:|:------------------------:|:---------------------------------------------:|
|Mud crab           | *Scylla paramamosain*    | C:72.0%[S:66.6%,D:5.4%],F:2.7%,M:25.3%,n:1013 |
|Black tiger shrimp | *Penaeus monodon*        | C:80.1%[S:77.5%,D:2.6%],F:5.5%,M:14.4%,n:1013 |
|Red clay crayfish  | *Cherax quadricarinatus* | C:83.7%[S:82.3%,D:1.4%],F:7.4%,M:8.9%,n:1013  |
|Brine shrimp       | *Artemia franciscana*    | C:81.0%[S:75.1%,D:5.9%],F:4.9%,M:14.1%,n:1013 |
*Note: additional species not available on RefSeq were also tested but excluded due to incompleteness*
  
  

