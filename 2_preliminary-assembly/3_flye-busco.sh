# Check completeness with BUSCO

module load busco/v5.3.1
busco -i all-runs_pacbio-raw/assembly.fasta -o flye_all10k-raw_busco -l /data/app/busco-5.3.1/lineages/arthropoda_odb10/ --mode genome --long -c 12
