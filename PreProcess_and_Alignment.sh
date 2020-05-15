#!/bin/bash

# Identify QC stats of all the fastqs with fastqc. 
for f1 in *.fastq.gz
do
	fastqc f1
done

# Check out the htmls generated for all the fastqc runs. 

# First, use UMI-tools to whitelist likely barcodes and extract them from the reads to move them onto the tag on the FASTQ's reads.
# Barcodes are on R2, so we perform whitelisting on R2s using the regex pattern shown. We then extract from R1. 

mkdir ExtractedFastqs
for f1 in *_R1_001.fastq.gz
do
	f2=${f1%%_R1_001.fastq.gz}"_R2_001.fastq.gz"
        umi_tools whitelist -I $f2 --bc-pattern="(?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT){s<=2}(?P<cell_2>.{8})(?P<umi_1>.{6})T{3}.*" --extract-method=regex --log2stderr "log_$f2" > whitelist.txt 
	umi_tools extract -I $f1 --bc-pattern2="(?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT){s<=2}(?P<cell_2>.{8})(?P<umi_1>.{6})T{3}.*" --extract-method=regex --read2-in="$f2" --stdout="extracted_$f1" --read2-out="extracted_$f2" --error-correct-cell --filter-cell-barcode --whitelist=whitelist.txt  
	mv "extracted_$f1" ExtractedFastqs
	mv "extracted_$f2" ExtractedFastqs
done

cd ExtractedFastqs

# Use STAR to align reads to the genome. The genome directory built by STAR needs to be in the current working directory.
# "GenomeDir" should be changed to the directory name that contains a STAR directory for a STAR-indexed GRCm38 genome.
# Genome can be downloaded from: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.p6.genome.fa.gz

for f1 in *_R1_001.fastq.gz
do
	f2=${f1%%_R1_001.fastq.gz}"_R2_001.fastq.gz"
	f3=${f1%%_R1_001.fastq.gz}
        STAR-2.7.3a/source/STAR --runThreadN 16 --genomeDir GenomeDir --readFilesIn $f1 $f2 --readFilesCommand zcat --outFileNamePrefix "$f3" --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate
done

# Use Subread's featureCounts to count transcripts mapped to each gene. 
# This command requires the GTF file for GRCm38 in directory. It can be downloaded at the link from which the genome was retrieved.
# Replace "GTF" with the GTF filename in the directory.

for f1 in *.bam
do
	subread-2.0.0-source/bin/featureCounts -a GTF -g gene_name -o "FC_$f1" -R BAM $f1 -T 16
	samtools sort "$f1.featureCounts.bam" -o "sorted_FC_$f1"
	samtools index "sorted_FC_$f1"
done

mkdir CountedFiles_Matrix

# Finally, we use UMI-tools to count the number of each transcript in each cell. Note that we count in "wide format".
# This outputs rows as genes and columns as cell names.

for f1 in *.bam
do
	f2=${f1%%.Aligned.sortedByCoord.out.bam}
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts -I $f1 -S "m_counts_$f2.tsv"
	mv "m_counts_$f2.tsv" CountedFiles_Matrix
done
