# GTEX_my_version
Building GTEX our version for Guis lab

First I wanted to get the docker image for GTEX. As docker is not supported on computecanada, I used
```
salloc --time=00:30:00 --mem=4G --account=def-ben
module load apptainer
apptainer pull gtex_rnaseq_V10.sif docker://broadinstitute/gtex_rnaseq:V10
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa gtex_refs/


get https://storage.cloud.google.com/gtex-resources/STAR_genomes/star_index_v19_oh75.tar.gz
