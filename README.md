# GTEX_my_version
Building GTEX our version for Guis lab

First I wanted to get the docker image for GTEX. As docker is not supported on computecanada, I used
```
salloc --time=00:30:00 --mem=4G --account=def-ben
module load apptainer
apptainer pull gtex_rnaseq_V10.sif docker://broadinstitute/gtex_rnaseq:V10
```
downloaded reference genome
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa gtex_refs/
```

copied the python script to remove ALTs and etc
```
with open('GRCh38.primary_assembly.genome.fa', 'r') as fasta:
    contigs = fasta.read()
contigs = contigs.split('>')
contig_ids = [i.split(' ', 1)[0] for i in contigs]

# exclude ALT, HLA and decoy contigs
filtered_fasta = '>'.join([c for i,c in zip(contig_ids, contigs)
    if not (i[-4:]=='_alt' or i[:3]=='HLA' or i[-6:]=='_decoy')])

with open('Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta', 'w') as fasta:
    fasta.write(filtered_fasta)
```
Then ran python vitual EN
```
module load python
virtualenv --no-download ~/ENV
source ~/ENV/bin/activate
pip install --no-index --upgrade pip
pip install numpy --no-index
```
and ran the script
```
python remove_ALT_HLA_Decoy.py
```
This creates "Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"



