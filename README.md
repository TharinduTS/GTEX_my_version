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
# These ateps are derived from GTEX github page "https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md"

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
Then downloaded ERCC spike-in reference annotations and unzipped them and processed as detailed below.
```
wget https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip
unzip ERCC92.zip
sed 's/ERCC-/ERCC_/g' ERCC92.fa > ERCC92.patched.fa
```
downloaded picard
```
wget https://github.com/broadinstitute/picard/releases/download/2.27.1/picard.jar
```

Then ERCC spike-in sequences were appended:
```
cat Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta ERCC92.patched.fa \
    > Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
```
The FASTA index and dictionary (required for Picard/GATK) were generated:
```
samtools faidx Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
java -jar picard.jar \
    CreateSequenceDictionary \
    R=Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta \
    O=Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict
```
# Reference annotation

downloaded annotation file ,gunzip and rename
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
gunzip gencode.v39.annotation.gtf.gz
mv gencode.v39.annotation.gtf gencode.v39.GRCh38.annotation.gtf
```
and collapse script
```
wget https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/gene_model/collapse_annotation.py
```
# installed bx-python
```
module load python/3.10
python -m venv gtex_env
source gtex_env/bin/activate
pip install --upgrade pip
pip install bx-python
```
For gene-level quantifications, the annotation was collapsed with the script used in the GTEx pipeline:
```
python collapse_annotation.py --stranded \
gencode.v39.GRCh38.annotation.gtf \
gencode.v39.GRCh38.genes.gtf

```

Gene- and transcript-level attributes were added to the ERCC GTF with the following Python code:
```add_attributes.py
with open('ERCC92.gtf') as exon_gtf, open('ERCC92.patched.gtf', 'w') as gene_gtf:
    for line in exon_gtf:
        f = line.strip().split('\t')
        f[0] = f[0].replace('-', '_')  # required for RNA-SeQC/GATK (no '-' in contig name)
        f[5] = '.'

        attr = f[8]
        if attr[-1] == ';':
            attr = attr[:-1]
        attr = dict([i.split(' ') for i in attr.replace('"', '').split('; ')])
        # add gene_name, gene_type
        attr['gene_name'] = attr['gene_id']
        attr['gene_type'] = 'ercc_control'
        attr['gene_status'] = 'KNOWN'
        attr['level'] = 2
        for k in ['id', 'type', 'name', 'status']:
            attr[f'transcript_{k}'] = attr[f'gene_{k}']

        attr_str = []
        for k in ['gene_id', 'transcript_id', 'gene_type', 'gene_status', 'gene_name',
            'transcript_type', 'transcript_status', 'transcript_name']:
            attr_str.append(f'{k} "{attr[k]}";')
        attr_str.append(f"level {attr['level']};")
        f[8] = ' '.join(attr_str)

        # write gene, transcript, exon
        gene_gtf.write('\t'.join(f[:2] + ['gene'] + f[3:]) + '\n')
        gene_gtf.write('\t'.join(f[:2] + ['transcript'] + f[3:]) + '\n')
        f[8] = ' '.join(attr_str[:2])
        gene_gtf.write('\t'.join(f[:2] + ['exon'] + f[3:]) + '\n')
```
run it 
```
python add_attributes.py
```
The ERCC annotation was appended to the reference GTFs:
```
cat gencode.v39.GRCh38.annotation.gtf ERCC92.patched.gtf \
    > gencode.v39.GRCh38.annotation.ERCC.gtf
cat gencode.v39.GRCh38.genes.gtf ERCC92.patched.gtf \
    > gencode.v39.GRCh38.ERCC.genes.gtf
```
had to patch ERCC with following
```patch_ERCC
with open("ERCC92.gtf") as exon_gtf, open("ERCC92.patched.gtf", "w") as gene_gtf:
    for line in exon_gtf:
        f = line.strip().split("\t")
        f[0] = f[0].replace("-", "_")  # STAR doesn't accept '-' in contig names
        f[5] = '.'
        attr = f[8].rstrip(";")
        attr = dict([i.split(" ") for i in attr.replace('"', '').split("; ")])
        attr["gene_name"] = attr["gene_id"]
        attr["gene_type"] = "ercc_control"
        attr["gene_status"] = "KNOWN"
        attr["level"] = "2"
        for k in ["id", "type", "name", "status"]:
            attr[f"transcript_{k}"] = attr[f"gene_{k}"]
        attr_str = [f'{k} "{attr[k]}";' for k in [
            "gene_id", "transcript_id", "gene_type", "gene_status",
            "gene_name", "transcript_type", "transcript_status", "transcript_name"
        ]]
        attr_str.append(f'level {attr["level"]};')
        f[8] = " ".join(attr_str)
        gene_gtf.write("\t".join(f[:2] + ["gene"] + f[3:]) + "\n")
        gene_gtf.write("\t".join(f[:2] + ["transcript"] + f[3:]) + "\n")
        f[8] = " ".join(attr_str[:2])
        gene_gtf.write("\t".join(f[:2] + ["exon"] + f[3:]) + "\n")
```
run it
```
python patch_ERCC
```
Append to GENCODE GTF
```
cat gencode.v39.GRCh38.annotation.gtf ERCC92.patched.gtf > gencode.v39.GRCh38.annotation.ERCC92.gtf
```

RSEM index
```
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load rsem
module load star

mkdir -p rsem_reference

rsem-prepare-reference \
--gtf gencode.v39.GRCh38.annotation.ERCC92.gtf \
--star \
--num-threads 4 \
Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta \
rsem_reference/rsem_reference
```

#*******Running the pipeline

I DOWNLOADED FASTQ FILES WITH SRA TOOLKIT. tHEREFORE SKIPPING THE BAM TO FASTQ CONVERSION BELOW
#*****
Comvert BAM to FASTQ
```BAM_to_FASTQ.sh
module load apptainer
apptainer run docker://broadinstitute/gtex_rnaseq     /bin/bash -c "rm -f read1_pipe read2_pipe read0_pipe && \
    /src/run_SamToFastq.py  /scratch/premacht/guis_lab/gtex_refs/samples/ENCFF069APB.bam -p ENCFF069APB -o ./samples"
```
#****

Download FASTQ with SRAtoolkit
```
module load sra-toolkit
fastq-dump --split-files --gzip SRR1553531
```
then test the read length so we can build star index accordingly
```
zcat samples/SRR1553531_1.fastq.gz | awk 'NR%4==2' | head -n 100 | awk '{print length($0)}' | sort | uniq -c
```
Mine was 101. Therefore I am building Star_index setting sjdbOverhang to 100 (it should be the length-1)

STAR index
The STAR index should be built to match the sequencing read length, specified by the sjdbOverhang parameter. These samples were sequenced using a 2x101 bp paired-end sequencing protocol, and the matching sjdbOverhang is 100.
CHANGE THE PLACES WITH 100 ACCORDINGLY      
```
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load star
mkdir -p star_index_oh100

STAR --runMode genomeGenerate \
--genomeDir star_index_oh100 \
--genomeFastaFiles Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta \
--sjdbGTFfile gencode.v39.GRCh38.annotation.ERCC92.gtf \
--sjdbOverhang 100 \
--runThreadN 4

```





