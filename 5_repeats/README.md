# Repeat annotation and masking

[RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) was used for de novo repeat identification and classification.  
To run RepeatModeler (v2.0.4), the following commands were used:  
```
cd JonahCrabGenome/Analyses/4_repeat-analysis
ln -s ../JC-genome_v3-230128.fasta .

module load repeatmodeler/v2.0.4
BuildDatabase -name jc-genome JC-genome_v3-230128.fasta
RepeatModeler -threads 32 -database jc-genome -LTRStruct
```
  
[RepeatMasker](https://www.repeatmasker.org/RepeatMasker/) then soft-masked the genome sequence using the repeat families identified with RepeatModeler.
```
module load repeatmasker/v4.1.4
RepeatMasker -pa 10 -lib ../jc-genome-families.fa -xsmall -gff JC-genome_v3-230128.fasta
```
  
Other programs used by RepeatModeler and RepeatMasker: NCBI/RMBlast (v2.13.0), TRF (v4.09.1), RECON (v1.08), RepeatScout (v1.0.6).  
Databases: Dfam (h5)
