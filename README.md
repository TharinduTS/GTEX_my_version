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
Their run_STAR.py script had bugs. So I had to fix it and save it as run_STAR_edited.py
```run_STAR_edited.py
#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import os
import subprocess
import gzip
import shutil
from datetime import datetime
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


parser = argparse.ArgumentParser(description='Run STAR')
parser.add_argument('index', help='Path to STAR index')
parser.add_argument('fastq', nargs='+', help='FASTQ input. Format: fastq1 [fastq2], or comma-separated lists for each if multiple FASTQs/mate.')
parser.add_argument('prefix', help='Prefix for output file names')
parser.add_argument('-o', '--output_dir', default='./', help='Output directory')
parser.add_argument('--annotation_gtf', default=None, help='Annotation in GTF format')
parser.add_argument('--outFilterMultimapNmax', default='20')
parser.add_argument('--alignSJoverhangMin', default='8')
parser.add_argument('--alignSJDBoverhangMin', default='1')
parser.add_argument('--outFilterMismatchNmax', default='999')
parser.add_argument('--outFilterMismatchNoverLmax', default='0.1')
parser.add_argument('--alignIntronMin', default='20')
parser.add_argument('--alignIntronMax', default='1000000')
parser.add_argument('--alignMatesGapMax', default='1000000')
parser.add_argument('--outFilterType', default='BySJout')
parser.add_argument('--outFilterScoreMinOverLread', default='0.33')
parser.add_argument('--outFilterMatchNmin', default='0')
parser.add_argument('--outFilterMatchNminOverLread', default='0.33')
parser.add_argument('--limitSjdbInsertNsj', default='1200000')
parser.add_argument('--outSAMstrandField', default='intronMotif')
parser.add_argument('--outFilterIntronMotifs', default='None', help="Use 'RemoveNoncanonical' for Cufflinks compatibility")
parser.add_argument('--alignSoftClipAtReferenceEnds', default='Yes')
parser.add_argument('--quantMode', default=['TranscriptomeSAM', 'GeneCounts'], nargs='+', help='Outputs read counts, and a BAM with reads in transcriptome coordinates')
parser.add_argument('--outSAMtype', default=['BAM', 'Unsorted'], nargs='+')
parser.add_argument('--outSAMunmapped', default='Within', help='Keep unmapped reads in output BAM')
parser.add_argument('--outSAMattrRGline', default=['ID:rg1', 'SM:sm1'], nargs='+', help='Adds read group line to BAM header; required by GATK')
parser.add_argument('--outSAMattributes', default=['NH', 'HI', 'AS', 'nM', 'NM', 'ch'], nargs='+')
parser.add_argument('--varVCFfile', default=None, help='VCF for the input sample; currently supports SNPs only')
parser.add_argument('--waspOutputMode', default='SAMtag')
parser.add_argument('--winAnchorMultimapNmax', default='50', help='int>0: max number of loci anchors are allowed to map to')
parser.add_argument('--genomeTransformOutput', default=None, nargs='+', help='string(s) which output to transform back to original genome')
parser.add_argument('--chimSegmentMin', default='15', help='Minimum chimeric segment length; switches on detection of chimeric (fusion) alignments')
parser.add_argument('--chimJunctionOverhangMin', default='15', help='Minimum overhang for a chimeric junction')
parser.add_argument('--chimOutType', default=['Junctions', 'WithinBAM', 'SoftClip'], nargs='+', help='')
parser.add_argument('--chimMainSegmentMultNmax', default='1', help='')
parser.add_argument('--chimOutJunctionFormat', default='0', help='Formatting for Chimeric.out.junction')
parser.add_argument('--genomeLoad', default='NoSharedMemory')
parser.add_argument('--sjdbFileChrStartEnd', default=None, help='SJ.out.tab file (e.g., from 1st pass). With this option, only one pass will be run')
parser.add_argument('--STARlong', action='store_true', help='Use STARlong instead of STAR')
parser.add_argument('-t', '--threads', default='4', help='Number of threads')
args = parser.parse_args()

if args.STARlong:
    starcmd = 'STARlong'
else:
    starcmd = 'STAR'

# set up command
cmd = f'{starcmd} --runMode alignReads --runThreadN {args.threads} --genomeDir {args.index}'
if args.annotation_gtf is not None:  # only needed if genome index was built w/o annotation
    cmd += f' --sjdbGTFfile {args.annotation_gtf}'
if args.sjdbFileChrStartEnd is None:
    cmd += ' --twopassMode Basic'
cmd += f' --outFilterMultimapNmax {args.outFilterMultimapNmax}'\
    + f' --alignSJoverhangMin {args.alignSJoverhangMin}'\
    + f' --alignSJDBoverhangMin {args.alignSJDBoverhangMin}'\
    + f' --outFilterMismatchNmax {args.outFilterMismatchNmax}'\
    + f' --outFilterMismatchNoverLmax {args.outFilterMismatchNoverLmax}'\
    + f' --alignIntronMin {args.alignIntronMin}'\
    + f' --alignIntronMax {args.alignIntronMax}'\
    + f' --alignMatesGapMax {args.alignMatesGapMax}'\
    + f' --outFilterType {args.outFilterType}'\
    + f' --outFilterScoreMinOverLread {args.outFilterScoreMinOverLread}'\
    + f' --outFilterMatchNmin {args.outFilterMatchNmin}'\
    + f' --outFilterMatchNminOverLread {args.outFilterMatchNminOverLread}'\
    + f' --limitSjdbInsertNsj {args.limitSjdbInsertNsj}'\
    + f" --readFilesIn {' '.join(args.fastq)}"
if args.fastq[0].endswith('.gz'):
    cmd += ' --readFilesCommand zcat'
cmd += f' --outFileNamePrefix {os.path.join(args.output_dir, args.prefix)}.'\
    + f' --outSAMstrandField {args.outSAMstrandField}'\
    + f' --outFilterIntronMotifs {args.outFilterIntronMotifs}'\
    + f' --alignSoftClipAtReferenceEnds {args.alignSoftClipAtReferenceEnds}'\
    + f" --quantMode {' '.join(args.quantMode)}"\
    + f" --outSAMtype {' '.join(args.outSAMtype)}"\
    + f' --outSAMunmapped {args.outSAMunmapped}'\
    + f' --genomeLoad {args.genomeLoad}'
if args.waspOutputMode=='SAMtag' and args.varVCFfile is not None:
    assert args.varVCFfile.endswith('.vcf.gz')
    cmd += f' --waspOutputMode SAMtag --varVCFfile <(zcat {args.varVCFfile})'
    if 'vW' not in args.outSAMattributes:
        args.outSAMattributes.append('vW')
        print("  * adding 'vW' tag to outSAMattributes", flush=True)
cmd += f' --winAnchorMultimapNmax {args.winAnchorMultimapNmax}'
if args.genomeTransformOutput is not None:
    cmd += f" --genomeTransformOutput {' '.join(args.genomeTransformOutput)}"
if int(args.chimSegmentMin)>0:
    cmd += f' --chimSegmentMin {args.chimSegmentMin}'\
        + f' --chimJunctionOverhangMin {args.chimJunctionOverhangMin}'\
        + f" --chimOutType {' '.join(args.chimOutType)}"\
        + f' --chimMainSegmentMultNmax {args.chimMainSegmentMultNmax}'\
        + f' --chimOutJunctionFormat {args.chimOutJunctionFormat}'
cmd += f" --outSAMattributes {' '.join(args.outSAMattributes)}"\
     + f" --outSAMattrRGline {' '.join(args.outSAMattrRGline)}"
if args.sjdbFileChrStartEnd is not None:
    cmd += f' --sjdbFileChrStartEnd {args.sjdbFileChrStartEnd}'

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# run STAR
subprocess.check_call(cmd, shell=True, executable='/bin/bash')

# postprocessing
with cd(args.output_dir):
    # set permissions
    for r,d,f in os.walk(f'{args.prefix}._STARpass1'):
        os.chmod(r, 0o755)

    # delete unneeded files
        # delete unneeded files
    shutil.rmtree(f'{args.prefix}._STARgenome')
    if os.path.exists(f'{args.prefix}._STARtmp'):
        shutil.rmtree(f'{args.prefix}._STARtmp')

    # sort BAM (use samtools to get around the memory gluttony of STAR)
    print(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] Sorting BAM', flush=True)
    cmd = f'samtools sort --threads {args.threads} -o {args.prefix}.Aligned.sortedByCoord.out.bam {args.prefix}.Aligned.out.bam'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    os.remove(f'{args.prefix}.Aligned.out.bam')
    print(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] Finished sorting BAM', flush=True)

    # index BAM
    print(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] Indexing BAM', flush=True)
    cmd = f'samtools index {args.prefix}.Aligned.sortedByCoord.out.bam'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    print(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] Finished indexing BAM', flush=True)

    # rename and compress outputs
    subprocess.check_call(f'gzip {args.prefix}.SJ.out.tab', shell=True, executable='/bin/bash')
    with cd(f'{args.prefix}._STARpass1'):
        os.rename('SJ.out.tab', f'{args.prefix}.SJ.pass1.out.tab')
        subprocess.check_call(f'gzip {args.prefix}.SJ.pass1.out.tab', shell=True, executable='/bin/bash')

    if os.path.exists(f'{args.prefix}.ReadsPerGene.out.tab'):
        subprocess.check_call(f'gzip {args.prefix}.ReadsPerGene.out.tab', shell=True, executable='/bin/bash')

    # sort and index chimeric BAM
    if os.path.exists(f'{args.prefix}.Chimeric.out.sam'):
        cmd = f'samtools sort --threads {args.threads} -o {args.prefix}.Chimeric.out.sorted.bam {args.prefix}.Chimeric.out.sam'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        cmd = f'samtools index {args.prefix}.Chimeric.out.sorted.bam'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        os.remove(f'{args.prefix}.Chimeric.out.sam')

    if os.path.exists(f'{args.prefix}.Chimeric.out.junction'):
        subprocess.check_call(f'gzip -f {args.prefix}.Chimeric.out.junction', shell=True, executable='/bin/bash')

```
Then ran it as following
First allocate some memory
```
salloc --cpus-per-task=4 --mem=64G --time=2:00:00
```
then run
```
sample_id=SRR1553531
apptainer run -B /scratch/premacht/guis_lab:/data /scratch/premacht/guis_lab/gtex_rnaseq_V10.sif     /bin/bash -c "python3 /data/run_STAR_edited.py \
        /data/gtex_refs/star_index_oh100 \
        /data/gtex_refs/samples/${sample_id}_1.fastq.gz \
        /data/gtex_refs/samples/${sample_id}_2.fastq.gz \
        ${sample_id} \
        --threads 4 \
        --output_dir /data/gtex_refs/samples/star_out_${sample_id}"
```
#*****THIS IS OPTIONAL. BUT I AM TRYING IT
getting picard and making mock BAM with fastqs
downloading the bam file to sync it
```
module load sra-toolkit
prefetch SRR1553531
```
convert to BAM
```
sam-dump SRR1553531 > SRR1553531.bam
```
then sync BAM
```
apptainer pull picard.sif docker://broadinstitute/picard
module load apptainer

sample_id=SRR1553531

apptainer run -B /scratch/premacht/guis_lab/gtex_refs:/data /scratch/premacht/guis_lab/gtex_rnaseq_V10.sif \
    /bin/bash -c "/src/run_bamsync.sh \
        /data/samples/${sample_id}.bam \
        /data/samples/star_out_${sample_id}.Aligned.sortedByCoord.out.patched.bam \
        /data/samples/star_out_${sample_id}"

```



this BAM header was missing a line
so had to add it
```
samtools view -H star_out_SRR1553531.Aligned.sortedByCoord.out.patched.bam > header.sam
echo -e "@RG\tID:SRR15.S\tSM:SRR1553531\tPL:ILLUMINA" >> header.sam
samtools reheader header.sam star_out_SRR1553531.Aligned.sortedByCoord.out.patched.bam > star_out_SRR1553531.Aligned.sortedByCoord.out.patched.fixed.bam
```
Mark duplicated
```
module load apptainer

sample_id=SRR1553531

apptainer run -B /scratch/premacht/guis_lab/gtex_refs:/data /scratch/premacht/guis_lab/gtex_rnaseq_V10.sif \
    /bin/bash -c "/src/run_MarkDuplicates.py \
        /data/samples/star_out_${sample_id}.Aligned.sortedByCoord.out.patched.fixed.bam \
        star_out_${sample_id}.Aligned.sortedByCoord.out.patched.md \
        --output_dir /data/samples"
```

RNA-SeQC v2 parses exon IDs to build its internal feature map. If the same exon ID appears more than once (even across different transcripts), it throws a fatal error.
So I had to collapse annotation
```
module load python/3.10
source gtex_env/bin/activate
pip install pandas

python collapse_annotation.py \
    --stranded \
    gencode.v39.GRCh38.annotation.gtf \
    gencode.v39.GRCh38.ERCC.genes.gtf
```

Then you can  RNA-SeQC





























sed 's/\(ENSE[0-9]*\)\.[0-9]\+/\1/g' gencode.v39.GRCh38.annotation.gtf > gencode.v39.GRCh38.annotation.noversion.gtf
```




testing
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


module load apptainer

apptainer run -B /scratch/premacht/guis_lab:/data /scratch/premacht/guis_lab/gtex_rnaseq_V10.sif \
    STAR \
    --runThreadN 4 \
    --genomeDir /data/gtex_refs/star_index_oh100 \
    --readFilesIn /data/gtex_refs/samples/SRR1553531_1.fastq.gz /data/gtex_refs/samples/SRR1553531_2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix /data/gtex_refs/samples/star_out_SRR1553531_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM
```

```
apptainer run -B /scratch/premacht/guis_lab:/data /scratch/premacht/guis_lab/gtex_rnaseq_V10.sif     /bin/bash -c "/src/run_bamsync.sh \
        ${input_bam} \
        /data/gtex_refs/samples/star_out_${sample_id}_Aligned.sortedByCoord.out.bam \
        /data/gtex_refs/samples/star_out_${sample_id}"
```
```
module load samtools

samtools index /scratch/premacht/guis_lab/gtex_refs/samples/star_out_SRR1553531_Aligned.sortedByCoord.out.bam
```













module load apptainer

sample_id=SRR1553531

apptainer run -B /scratch/premacht/guis_lab:/data /scratch/premacht/guis_lab/gtex_rnaseq_V10.sif \
    /bin/bash -c "/src/run_STAR.py \
        /data/gtex_refs/star_index_oh100 \
        /data/gtex_refs/samples/${sample_id}_1.fastq.gz \
        /data/gtex_refs/samples/${sample_id}_2.fastq.gz \
        ${sample_id} \
        --threads 4 \
        --output_dir /tmp/star_out && mv /tmp/star_out /data/gtex_refs/samples/star_out_${sample_id}"

```

