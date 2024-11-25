#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status
# Uncomment the following line to print commands and their arguments as they are executed
# set -x

# Set up Spack
source /n/home00/abakirbas/Desktop/spack/share/spack/setup-env.sh || { echo "Failed to source Spack"; exit 1; }

# Create and activate Spack environment (if it doesn't exist)
if ! spack env list | grep -q rnaseq_env; then
    spack env create rnaseq_env || { echo "Failed to create Spack environment"; exit 1; }
fi
eval "$(spack env activate --sh rnaseq_env)" || { echo "Failed to activate Spack environment"; exit 1; }

# Install packages (if not already installed)
if ! spack find fastqc hisat2 stringtie samtools transdecoder > /dev/null 2>&1; then
    spack add fastqc hisat2 stringtie samtools transdecoder || { echo "Failed to add packages"; exit 1; }
    spack install || { echo "Failed to install packages"; exit 1; }
fi

# Load packages
spack load fastqc hisat2 stringtie samtools transdecoder || { echo "Failed to load packages"; exit 1; }

# Set variables
GENOME="/n/netscratch/davis_lab/Everyone/abakirbas/b_caapi.p_ctg.fa"
READS_DIR="/n/netscratch/davis_lab/Everyone/abakirbas/Caapi_RNA_20221216"
OUTPUT_DIR="/n/netscratch/davis_lab/Everyone/abakirbas"
SCRIPTS_DIR="/n/netscratch/davis_lab/Everyone/abakirbas/scripts"

# Create output and scripts directories
mkdir -p "$OUTPUT_DIR" "$SCRIPTS_DIR"
mkdir -p "${OUTPUT_DIR}/fastqc"
mkdir -p "${OUTPUT_DIR}/aligned"
mkdir -p "${OUTPUT_DIR}/assembled"
mkdir -p "${OUTPUT_DIR}/logs"

# Create a sample list with full file names
find "${READS_DIR}" -name "*_R1_001.fastq.gz" | sort > "${SCRIPTS_DIR}/sample_list_R1.txt"
find "${READS_DIR}" -name "*_R2_001.fastq.gz" | sort > "${SCRIPTS_DIR}/sample_list_R2.txt"
paste "${SCRIPTS_DIR}/sample_list_R1.txt" "${SCRIPTS_DIR}/sample_list_R2.txt" > "${SCRIPTS_DIR}/sample_list_full.txt"

# Clear previous job IDs file
> "${SCRIPTS_DIR}/sample_job_ids.txt"

# Generate individual sample job scripts
while read R1 R2; do
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    cat << EOF > "${SCRIPTS_DIR}/process_${SAMPLE}.sh"
#!/bin/bash
#SBATCH --job-name=rnaseq_${SAMPLE}
#SBATCH --output="${OUTPUT_DIR}/logs/rnaseq_${SAMPLE}_%j.out"
#SBATCH --error="${OUTPUT_DIR}/logs/rnaseq_${SAMPLE}_%j.err"
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --partition=sapphire

# Load packages
source /n/home00/abakirbas/Desktop/spack/share/spack/setup-env.sh
# Activate Spack environment in the current shell
eval "\$(spack env activate --sh rnaseq_env)"
spack load fastqc hisat2 stringtie samtools

# Set variables
GENOME="${GENOME}"
OUTPUT_DIR="${OUTPUT_DIR}"
R1_FILE="${R1}"
R2_FILE="${R2}"
SAMPLE="${SAMPLE}"

# Verify SLURM_CPUS_PER_TASK
echo "SLURM_CPUS_PER_TASK = \$SLURM_CPUS_PER_TASK"
if [ -z "\$SLURM_CPUS_PER_TASK" ]; then
    SLURM_CPUS_PER_TASK=20  # Set to the value specified in SBATCH directives
    echo "SLURM_CPUS_PER_TASK was unset. Setting it to \$SLURM_CPUS_PER_TASK."
fi

# Check if HISAT2 index exists, if not create it
if [ ! -f "\${GENOME%.*}_index.1.ht2" ]; then
    echo "Creating HISAT2 index..."
    hisat2-build "\$GENOME" "\${GENOME%.*}_index"
fi

# FastQC
fastqc "\$R1_FILE" "\$R2_FILE" -o "\${OUTPUT_DIR}/fastqc"

# HISAT2 alignment
hisat2 -p "\$SLURM_CPUS_PER_TASK" \
    -x "\${GENOME%.*}_index" \
    -1 "\$R1_FILE" \
    -2 "\$R2_FILE" \
    -S "\${OUTPUT_DIR}/aligned/\${SAMPLE}.sam"

# Check if SAM file was created
if [ ! -f "\${OUTPUT_DIR}/aligned/\${SAMPLE}.sam" ]; then
    echo "Error: SAM file was not created. Exiting."
    exit 1
fi

# Convert SAM to BAM, sort, and index
samtools view -bS "\${OUTPUT_DIR}/aligned/\${SAMPLE}.sam" | \
    samtools sort -o "\${OUTPUT_DIR}/aligned/\${SAMPLE}.sorted.bam"
samtools index "\${OUTPUT_DIR}/aligned/\${SAMPLE}.sorted.bam"

# StringTie assembly and quantification
stringtie "\${OUTPUT_DIR}/aligned/\${SAMPLE}.sorted.bam" \
    -o "\${OUTPUT_DIR}/assembled/\${SAMPLE}.gtf" \
    -p "\$SLURM_CPUS_PER_TASK"

# Remove intermediate SAM file to save space
rm "\${OUTPUT_DIR}/aligned/\${SAMPLE}.sam"

# Deactivate Spack environment
eval "\$(spack env deactivate --sh)"

echo "Processing of \${SAMPLE} complete."
EOF

    # Submit the job and store the job ID
    JOB_ID=$(sbatch "${SCRIPTS_DIR}/process_${SAMPLE}.sh" | awk '{print $4}')
    if [ $? -eq 0 ]; then
      echo "$JOB_ID" >> "${SCRIPTS_DIR}/sample_job_ids.txt"
      echo "Submitted job $JOB_ID for sample ${SAMPLE}"
    else
      echo "Failed to submit job for sample ${SAMPLE}"
    fi
done < "${SCRIPTS_DIR}/sample_list_full.txt"

# Deactivate Spack environment
eval "$(spack env deactivate --sh)"