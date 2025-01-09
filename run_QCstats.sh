#!/bin/bash

# Remove any existing QCstats.txt file to start fresh
rm -f QCstats.txt

# A comment line and date/time for time stamping
echo "# QC Statistics generated on $(date)" >> QCstats.txt
echo "# This file contains counts of various features from annotation files" >> QCstats.txt
echo "" >> QCstats.txt

# Counting the number of genes from the TransDecoder genome-based GFF3 file
echo "# Number of genes predicted in bcaapi_merged_transcripts.fa.transdecoder.genome.gff3:" >> QCstats.txt
grep -c $'\tgene\t' bcaapi_merged_transcripts.fa.transdecoder.genome.gff3 >> QCstats.txt
echo "" >> QCstats.txt

# Counting the number of genes from the integrated Transdecoder annotation file
echo "# Number of genes predicted (after merging ncRNAs with Transdecoder output) in bcaapi_merged_gffread.gff3_wCDSfeatures.gff3:" >> QCstats.txt
grep -c $'\tgene\t' bcaapi_merged_gffread.gff3_wCDSfeatures.gff3 >> QCstats.txt
echo "" >> QCstats.txt

# Counting the number of predicted ORFs from the transdecoder cds file
echo "# Number of predicted ORFs in bcaapi_merged_transcripts.fa.transdecoder.cds:" >> QCstats.txt
grep -c "^>" bcaapi_merged_transcripts.fa.transdecoder.cds >> QCstats.txt
echo "" >> QCstats.txt

# Counting the number of transcripts from the raw StringTie GTF file
echo "# Number of transcripts from bcaapi_merged.gtf:" >> QCstats.txt
grep -c $'\ttranscript\t' bcaapi_merged.gtf >> QCstats.txt
echo "" >> QCstats.txt

# Count the number of predicted proteins from the TransDecoder peptide file
echo "# Number of predicted proteins in bcaapi_merged_transcripts.fa.transdecoder.pep:" >> QCstats.txt
grep -c "^>" bcaapi_merged_transcripts.fa.transdecoder.pep >> QCstats.txt

echo "" >> QCstats.txt
echo "# End of QC Statistics" >> QCstats.txt