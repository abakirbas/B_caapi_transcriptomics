# Load required library
library(Biostrings)

# setwd("~/Desktop/Davis Lab/B_caapi_transcriptomics/")
fasta_file <- "./bcaapi_merged_transcripts.fa.transdecoder.pep"

# Read nucleotide sequences
# seqs <- readDNAStringSet(fasta_file)

# For proteins, use readAAStringSet
seqs <- readAAStringSet(fasta_file)

# Extract lengths
seq_lengths <- width(seqs)

# Combine names and lengths, but the names extracted by readAAStringSet are not clean
# I looked at the file and determined a pattern to extract the substring that can be used as the "seq_name"
# I'm using "look behind" and "look ahead" to cleanly extract the substring e.g., "MSTRG.1".
# (?<=GENE\\.) means match after 'GENE.'
# [^~]* matches all non-tilde chars
# (?=~~) means match before '~~'
result <- data.frame(
  seq_name = stringr::str_extract(names(seqs), "(?<=GENE\\.)[^~]*(?=~~)"),
  length = seq_lengths
)

# Write results to a tab-delimited file
write.table(result, file = paste0(fasta_file, ".lengths.txt"), sep="\t", quote=FALSE, row.names=FALSE)

# Read results from file
result <- read.table("./bcaapi_merged_transcripts.fa.transdecoder.pep.lengths.txt", header=TRUE, sep="\t")

# Calculate the bin width using Scott's rule
bin_width_scott <- round(3.49 * sd(result$length) / length(result$length)^(1/3))

# Plot the distribution as a histogram using ggpubr package
library(ggpubr)

gghistogram(result$length, bins = bin_width_scott, xlab = "Protein Length (log10 scale)", ylab = "Frequency", fill = "darkred", alpha = 0.3
            ) + 
  scale_x_log10() +
theme_bw() +
  ggtitle("B. caapi protein Length Distribution") +
  #remove background
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size = 16)) + 
  # Add mean value
  geom_vline(aes(xintercept = mean(result$length)), color = "blue", linetype = "dashed", size = 1)

# Create histogram
ggplot(result, aes(x=length)) +
  geom_histogram(bins=50, color="black", fill="blue") +
  theme_bw() +
  labs(title="Protein Length Distribution", x="Length (Amino Acids)", y="Count")

hist_file <- paste0(fasta_file, ".length_histogram.png")
ggsave(hist_file, plot=p, width=8, height=6, dpi=300)