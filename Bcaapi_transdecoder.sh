#!/bin/bash
#SBATCH --job-name=tdecoder
#SBATCH --output="/n/netscratch/davis_lab/Everyone/abakirbas/logs/tdecoder_%j.out"
#SBATCH --error="/n/netscratch/davis_lab/Everyone/abakirbas/logs/tdecoder_%j.err"
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --partition=sapphire

# Load Mamba/Conda
module purge
module load python
source /n/sw/Miniforge3-24.7.1-0/etc/profile.d/mamba.sh

# Activate the rnaseq environment (where stringtie, transdecoder and blast are installed)
conda init
conda activate rnaseq

# Set variables
OUTPUT_DIR="/n/netscratch/davis_lab/Everyone/abakirbas"
ASSEMBLED_DIR="${OUTPUT_DIR}/assembled"
MERGED_GTF="${OUTPUT_DIR}/bcaapi_merged.gtf"
GENOME_FA="${OUTPUT_DIR}/b_caapi.p_ctg.fa"
TRANSCRIPT_FA="${OUTPUT_DIR}/bcaapi_merged_transcripts.fa"
TRANSCRIPTS_GFF3="${OUTPUT_DIR}/bcaapi_merged.gff3"
BLASTP_HITS="${OUTPUT_DIR}/blastp_hits.outfmt6"
PROTEIN_DB_FA="${OUTPUT_DIR}/uniprot_sprot.fasta"  # Path to the protein database FASTA file
BLASTP_DB="${OUTPUT_DIR}/protein_db"  # Name of the BLAST database to be created

# Verify SLURM_CPUS_PER_TASK
echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
if [ -z "$SLURM_CPUS_PER_TASK" ]; then
    SLURM_CPUS_PER_TASK=40  # Set to the value specified in SBATCH directives
    echo "SLURM_CPUS_PER_TASK was unset. Setting it to $SLURM_CPUS_PER_TASK."
fi

# Merge GTF files
echo "Merging GTF files with StringTie..."
stringtie --merge -p "$SLURM_CPUS_PER_TASK" \
    -o "$MERGED_GTF" \
    "${ASSEMBLED_DIR}"/*.gtf

cd "$OUTPUT_DIR"

# Step 1: Convert GTF to transcript FASTA
echo "Converting GTF to transcript FASTA..."
gtf_genome_to_cdna_fasta.pl "$MERGED_GTF" "$GENOME_FA" > "$TRANSCRIPT_FA"

# Step 2: Convert GTF to alignment GFF3
echo "Converting GTF to alignment GFF3..."
gtf_to_alignment_gff3.pl "$MERGED_GTF" > "$TRANSCRIPTS_GFF3"

# Step 3: Run TransDecoder.LongOrfs
echo "Running TransDecoder.LongOrfs..."
TransDecoder.LongOrfs -t "$TRANSCRIPT_FA"

# Step 4: Prepare BLASTP database
echo "Creating BLASTP database..."
makeblastdb -in "$PROTEIN_DB_FA" -dbtype prot -out "$BLASTP_DB"

# Step 5: Optional BLASTP search for homology evidence, evalue set to 1e-4 to minimize the filtering out of real orfs per Freedman and Sackton (2024)
echo "Running BLASTP search..."
blastp -query "${TRANSCRIPT_FA}.transdecoder_dir/longest_orfs.pep" \
       -db "$BLASTP_DB" \
       -outfmt 6 \
       -evalue 1e-4 \
       -num_threads "$SLURM_CPUS_PER_TASK" \
       > "$BLASTP_HITS"

# Step 6: Run TransDecoder.Predict with BLASTP hits
echo "Running TransDecoder.Predict..."
TransDecoder.Predict -t "$TRANSCRIPT_FA" \
    --single_best_only \
    --retain_blastp_hits "$BLASTP_HITS"

# Step 7: Map ORFs back to genome coordinates
echo "Mapping ORFs back to genome coordinates..."
cdna_alignment_orf_to_genome_orf.pl "${TRANSCRIPT_FA}.transdecoder.gff3" \
    "$TRANSCRIPTS_GFF3" \
    "$TRANSCRIPT_FA" \
    > "${TRANSCRIPT_FA}.transdecoder.genome.gff3"

echo "TransDecoder job is complete!"

# Deactivate the rnaseq environment
conda deactivate