!/usr/bin/env bash


#Downloading the GRCh38 unpatched genome and assembly report
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
/home/user/bowtie2-2.2.2/bowtie2-build $PWD/data/GRCh38_alignment.fa $PWD/data/GRCh38_index

#Getting .bam files for analysis
for i in {6792585..6792587}; do

    #Downloading the ChIP-Seq data
    prefetch SRR${i}
    fastq-dump SRR${i} -O $PWD/data

    #Align ChIP-Seq reads to the GRCh38 reference genome
    /home/user/bowtie2-2.2.2/bowtie2 -x $PWD/data/GRCh38_index -U $PWD/data/SRR${i}.fastq -S $PWD/data/SRR${i}.sam

    #Converting .sam to .bam and sorting alignments by name
    samtools sort $PWD/data/SRR${i}.sam -o $PWD/data/bams/SRR${i}_sorted.bam

    #Removing PCR duplicates
    samtools markdup -r $PWD/data/bams/SRR${i}_sorted.bam $PWD/data/bams/SRR${i}_markdup.bam

    #Indexing the .bam files
    samtools index $PWD/data/bams/SRR${i}_markdup.bam
done

#Remove sorted .bam files
rm $PWD/data/bams/*_sorted.bam

#Prompt user input
while true; do
    read -p "Please provide the SRR Control data ID ("None" if not available): " ctrl
    case $ctrl in
        [SRR]* ) mv $PWD/data/bams/${ctrl}_markdup.bam* $PWD/data; exit;;
        None ) exit;;
        * ) echo "Please answer with SRR* or None";;
    esac
done
