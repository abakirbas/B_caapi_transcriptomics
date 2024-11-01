#!/bin/bash

# Set variables
GENOME="/path/to/your/denovo_assembled_genome.fasta"
READS_DIR="/path/to/your/reads_directory"
OUTPUT_DIR="/path/to/your/output_directory"
SCRIPTS_DIR="/path/to/your/scripts_directory"

# Create output and scripts directories
mkdir -p $OUTPUT_DIR $SCRIPTS_DIR

# Create a sample list
ls ${READS_DIR}/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' > ${SCRIPTS_DIR}/sample_list.txt

# Generate individual sample job scripts
while read sample; do
    cat << EOF > ${SCRIPTS_DIR}/process_${sample##*/}.sh
#!/bin/bash
#SBATCH --job-name=rnaseq_${sample##*/}
#SBATCH --output=${OUTPUT_DIR}/logs/rnaseq_${sample##*/}_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/rnaseq_${sample##*/}_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=your_partition_name

# Load necessary modules (adjust these based on your HPC environment)
module load fastqc
module load hisat2
module load stringtie
module load samtools

# Set variables
GENOME="$GENOME"
READS_DIR="$READS_DIR"
OUTPUT_DIR="$OUTPUT_DIR"
SAMPLE="${sample##*/}"

# FastQC
fastqc ${READS_DIR}/${SAMPLE}_R1.fastq.gz ${READS_DIR}/${SAMPLE}_R2.fastq.gz -o ${OUTPUT_DIR}/fastqc

# HISAT2 alignment
hisat2 -p \$SLURM_CPUS_PER_TASK \
    -x \${GENOME%.*}_index \
    -1 \${READS_DIR}/\${SAMPLE}_R1.fastq.gz \
    -2 \${READS_DIR}/\${SAMPLE}_R2.fastq.gz \
    -S \${OUTPUT_DIR}/\${SAMPLE}.sam

# Convert SAM to BAM, sort, and index
samtools view -bS \${OUTPUT_DIR}/\${SAMPLE}.sam | \
    samtools sort -o \${OUTPUT_DIR}/\${SAMPLE}.sorted.bam
samtools index \${OUTPUT_DIR}/\${SAMPLE}.sorted.bam

# StringTie assembly and quantification
stringtie \${OUTPUT_DIR}/\${SAMPLE}.sorted.bam \
    -o \${OUTPUT_DIR}/\${SAMPLE}.gtf \
    -p \$SLURM_CPUS_PER_TASK

# Remove intermediate SAM file to save space
rm \${OUTPUT_DIR}/\${SAMPLE}.sam

EOF

    # Submit the job and store the job ID
    JOB_ID=$(sbatch ${SCRIPTS_DIR}/process_${sample##*/}.sh | awk '{print $4}')
    echo $JOB_ID >> ${SCRIPTS_DIR}/sample_job_ids.txt
done < ${SCRIPTS_DIR}/sample_list.txt

# Create and submit the merge and TransDecoder job
cat << EOF > ${SCRIPTS_DIR}/merge_and_transdecoder.sh
#!/bin/bash
#SBATCH --job-name=merge_transdecoder
#SBATCH --output=${OUTPUT_DIR}/logs/merge_transdecoder_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/merge_transdecoder_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=your_partition_name

# Load necessary modules
module load stringtie
module load transdecoder

# Merge GTF files
stringtie --merge -p \$SLURM_CPUS_PER_TASK \
    -o ${OUTPUT_DIR}/merged.gtf \
    ${OUTPUT_DIR}/*.gtf

# TransDecoder to identify coding regions
cd $OUTPUT_DIR
TransDecoder.LongOrfs -t merged.gtf
TransDecoder.Predict -t merged.gtf

echo "RNA-seq analysis complete!"
EOF

# Submit the merge and TransDecoder job with dependency on all sample jobs
DEPEND_STRING=$(tr '\n' ':' < ${SCRIPTS_DIR}/sample_job_ids.txt | sed 's/:$//')
sbatch --dependency=afterok:$DEPEND_STRING ${SCRIPTS_DIR}/merge_and_transdecoder.sh

# Index the genome (only needs to be done once)
if [ ! -f ${GENOME}.1.ht2 ]; then
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --output=${OUTPUT_DIR}/logs/hisat2_index_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/hisat2_index_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=your_partition_name

module load hisat2

hisat2-build $GENOME ${GENOME%.*}_index
EOF
fi