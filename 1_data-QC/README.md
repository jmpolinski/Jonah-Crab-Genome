# QC of PacBio data for genome assembly

Data stats prior to QC:
```
module load seqkit
seqkit stats ./JC504_DTG-PacBio_all_raw.fastq
```
Output:
```
file                            format  type    num_seqs          sum_len  min_len   avg_len  max_len
JC504_DTG-PacBio_all_raw.fastq  FASTQ   DNA   22,144,746  345,729,859,605       50  15,612.3  285,394

```

Preliminary attempts to run some assemblers resulted in errors due to insufficient memory. Because of this, Trimmomatic was used to remove sequences less than 10 kb.

```
java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 16 -phred33 JC504_DTG-PacBio_all_raw.fastq JC504_DTG-PacBio_all_10kb.fastq MINLEN:10000
seqkit stats JC504_DTG-PacBio_all_10kb.fastq
``` 
Output
```
file                            format  type    num_seqs          sum_len  min_len   avg_len  max_len
JC504_DTG-PacBio_all_10kb.fastq  FASTQ   DNA   13,090,793  307,146,803,029   10,000  23,462.8  285,394
```
