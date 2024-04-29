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
  
  

## OrthoFinder ortholog identification
[OrthoFinder](https://github.com/davidemms/OrthoFinder) was used to identify groups of orthologous genes ("orthogroups") and generate a preliminary phylogenetic tree.   
See ```2_orthofinder.sh```  

## Time-scaled phylogenetic tree 
I used the [MEGA Timetree Wizard](https://www.megasoftware.net/web_help_11/Part_I_Getting_Started/A_Walk_Through_MEGA/Constructing_a_Timetree_(ML).htm) to make the phylogeneteic tree time-scaled.
To do this, I needed a tree, a sequence alignment, and time calibration points.  

Protein alignment: I used MUSCLE v5.1 to generate protein alignments for all single-copy orthologs identified by OrthoFinder, and then used [Gblocks](https://home.cc.umanitoba.ca/~psgendb/doc/Castresana/Gblocks_documentation.html) to extract conserved, well-aligned blocks from each of these and put them in a concatenated file. 
Code can be found in ```3_time-tree-alignments.sh```.   

Tree: I used the OrthoFinder tree with no branch lengths or confidence values (below)
```
((Bombus,Dmel),((Dpulex,Dmagna),((Lsalmonis,Tcalifornicus),(Hazteca,((Esinensis,(Ptrituberculatus,Cborealis)),((Hamericanus,Pclarkii),(Pjaponicus,(Lvannamei,Pchinensis))))))));
```
In MEGA (Windows GUI):
- Select Clock --> Compute Time Tree --> RelTime-ML
- The alignment (saved as .meg file) from above was input as the sequence data
- The tree above was input as the tree file
- The branch with D. melanogaster and B. pyrosoma was selected as the outgroup
- Calibration nodes were as follows: C. borealis - P. trituberculatus (max 225.0 MYA), L. vannamei - P. chinensis (>57.8<108.3; uniform distribution), H. americanus - P. clarkii (>241.0<321.6; uniform dist.), D. magna - T. californicus (275-541, uniform dist.) (from [TimeTree.org](timetree.org))
- Analysis preferences all default (JTT model; uniform rates among sites; use all sites)

## CAFE gene family expansion/contraction
The OrthoFinder output count table (in Results/Orthogroups/Orthogroups.GeneCount.tsv) was reformatted for CAFE:
- renamed first column "Desc"
- added new second column "Family ID", which was just the number part of the Orthogroup ID (i.e. OG0000001 = 1; note that OG0000000 could not be 0 so was highest number)
- saved as "counts4cafe_2024-04.txt"  
Checked if dataset needed filtering:
```
module load cafe/v4.2.1
python /data/resources/app_modules/CAFE/python_scripts/clade_size_filter.py -i counts4cafe_2024-04.txt -o filtered_cafe_input_2024-04.txt -s
# "No filtering was done!" output
```
Modified CAFE scipt and ran: ```./cafe_20240429.sh```
