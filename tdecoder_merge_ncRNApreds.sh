#!/bin/bash
#SBATCH --job-name=tdecoder_debug
#SBATCH --output="/n/netscratch/davis_lab/Everyone/abakirbas/logs/tdecoder_debug_%j.out"
#SBATCH --error="/n/netscratch/davis_lab/Everyone/abakirbas/logs/tdecoder_debug_%j.err"
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --partition=sapphire

# Load Mamba/Conda
module purge
module load python
source /n/sw/Miniforge3-24.7.1-0/etc/profile.d/mamba.sh

# Activate the rnaseq environment (where stringtie, transdecoder and gffread are installed)
mamba init
mamba activate rnaseq

# Set variables for input files
TDEC_GFF3="/n/netscratch/davis_lab/Everyone/abakirbas/bcaapi_merged_transcripts.fa.transdecoder.genome.gff3"
STRINGTIE_GTF="/n/netscratch/davis_lab/Everyone/abakirbas/bcaapi_merged.gtf"

# Run the Python script with the appropriate arguments
python MergeNcRnaPredstoProteinCodingPreds.py -tdecgff3 "$TDEC_GFF3" -gtf "$STRINGTIE_GTF"

# Deactivate the rnaseq environment after the job completes
conda deactivate

python MergeNcRnaPredstoProteinCodingPreds.py -tdecgff3 bcaapi_merged_transcripts.fa.transdecoder.genome.gff3 -gtf bcaapi_merged.gtf