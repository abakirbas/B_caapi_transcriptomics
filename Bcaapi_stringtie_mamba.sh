#!/bin/bash
#SBATCH --job-name=stringtie_analysis
#SBATCH --output=/n/netscratch/davis_lab/Everyone/abakirbas/logs/stringtie_%j.out
#SBATCH --error=/n/netscratch/davis_lab/Everyone/abakirbas/logs/stringtie_%j.err
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --partition=sapphire

# Load Mamba/Conda
module purge
module load python
source /n/sw/Miniforge3-24.7.1-0/etc/profile.d/mamba.sh

# Activate the rnaseq environment
conda init
conda activate rnaseq

# Set variables
OUTPUT_DIR="/n/netscratch/davis_lab/Everyone/abakirbas"
ALIGNED_DIR="${OUTPUT_DIR}/aligned"
ASSEMBLED_DIR="${OUTPUT_DIR}/assembled"
SCRIPTS_DIR="/n/netscratch/davis_lab/Everyone/abakirbas/scripts"
SAMPLE_LIST="${ALIGNED_DIR}/sample_list.txt"

# Create assembled directory if it doesn't exist
mkdir -p "${ASSEMBLED_DIR}"

# Verify SLURM_CPUS_PER_TASK
echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
if [ -z "$SLURM_CPUS_PER_TASK" ]; then
    SLURM_CPUS_PER_TASK=20  # Set to a default value if unset
    echo "SLURM_CPUS_PER_TASK was unset. Setting it to $SLURM_CPUS_PER_TASK."
fi

# Process each sample
while read SAMPLE; do
    echo "Processing sample: ${SAMPLE}"
    
    # Check if input BAM file exists
    if [ ! -f "${ALIGNED_DIR}/${SAMPLE}.sorted.bam" ]; then
        echo "Error: BAM file for ${SAMPLE} not found. Skipping."
        continue
    fi

    # StringTie assembly and quantification
    stringtie "${ALIGNED_DIR}/${SAMPLE}.sorted.bam" \
        -o "${ASSEMBLED_DIR}/${SAMPLE}.gtf" \
        -p "$SLURM_CPUS_PER_TASK"

    if [ $? -eq 0 ]; then
        echo "StringTie analysis completed for ${SAMPLE}"
    else
        echo "Error: StringTie analysis failed for ${SAMPLE}"
    fi
done < "$SAMPLE_LIST"

# Deactivate the rnaseq environment
conda deactivate

echo "StringTie analysis complete for all samples."