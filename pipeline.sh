#!/usr/bin/env bash


#Download the GRCh38 unpatched genome and assembly report
wget -P $PWD/data https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
wget -P $PWD/data https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_assembly_report.txt

#Unzipping the GRCh38 genome
gunzip $PWD/data/GCF_000001405.26_GRCh38_genomic.fna.gz

#Select entries in the genome FASTA file, belonging to the Primary Assembly or Mitochondrial chromosome
sort -k1,1V $PWD/data/GCF_000001405.26_GRCh38_assembly_report.txt | awk -F "\t" '$8 == "Primary Assembly" || $8 == "non-nuclear" {print $7}' > $PWD/data/subset_ids.txt
samtools faidx $PWD/data/GCF_000001405.26_GRCh38_genomic.fna -r $PWD/data/subset_ids.txt -o $PWD/data/GRCh38_subset.fa

#Replacing the RefSeq-style fasta headers with UCSC-style headers
awk -v FS="\t" 'NR==FNR {header[">"$7] = ">"$10; next} $0 ~ "^>" {sub($0, header[$0]); print}1' $PWD/data/GCF_000001405.26_GRCh38_assembly_report.txt $PWD/data/GRCh38_subset.fa > $PWD/data/GRCh38_alignment.fa

#Creating an index for GRCh38 reference genome
bowtie2-build $PWD/data/GRCh38_alignment.fa $PWD/data/GRCh38_index