from pathlib import Path


configfile: "config/example.yaml"


# Get data from configuration file
DATA_DIR = config["dir"]["data_dir"]
OUT_DIR = config["dir"]["output_dir"]

# Initialising variables
GENOME_FILE = Path(config["data"]["genome"]).name
ASSEMBLY_FILE = Path(config["data"]["assembly_report"]).name
CONTROLS_IDS = config["data"]["chipseq"]["control_ids"]
TARGET_IDS = config["data"]["chipseq"]["target_ids"]
ALL_IDS = CONTROLS_IDS + TARGET_IDS

# MACS2 parameters
BROAD_CUTOFF = config["parameters"]["macs2"]["broad_cutoff"]
GENOME_SIZE = config["parameters"]["macs2"]["genome_size"]
EXTENSION_SIZE = config["parameters"]["macs2"]["extension_size"]

# Define wildcards
index_suffixes = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]


rule all:
    """
    Defining the final expected output files.
    """
    input:
        f"{DATA_DIR}/{GENOME_FILE[:-3]}",
        f"{DATA_DIR}/{ASSEMBLY_FILE}", 
        f"{DATA_DIR}/genome_alignment.fa",
        [f"{DATA_DIR}/genome_index.{suffix}" for suffix in index_suffixes],
        expand(f"{DATA_DIR}/{{srr}}.fastq", srr=ALL_IDS),
        expand(f"{DATA_DIR}/bams/{{srr}}_markdup.bam", srr=CONTROLS_IDS),
        expand(f"{DATA_DIR}/bams/{{srr}}_markdup.bam.bai", srr=CONTROLS_IDS),
        expand(f"{DATA_DIR}/{{srr}}_markdup.bam", srr=TARGET_IDS) +
        expand(f"{DATA_DIR}/{{srr}}_markdup.bam.bai", srr=TARGET_IDS),
        expand(f"{DATA_DIR}/macs2/{{srr}}_peaks.broadPeak", srr=TARGET_IDS),
        expand(f"{DATA_DIR}/{{srr}}_constituent_enhancers.gff3", srr=TARGET_IDS),

# ----------------------------------- #
# 01. Reference Genome Preparation    #
# ----------------------------------- #

rule download_genome:
    """
    Download genome and assembly report.
    """
    output:
        genome = temp(f"{DATA_DIR}/{GENOME_FILE}"),
        assembly_report = f"{DATA_DIR}/{ASSEMBLY_FILE}"
    params:
        genome_url = config["data"]["genome"],
        assembly_report_url = config["data"]["assembly_report"],
    log:
        f"{OUT_DIR}/log/01_download_genome.log"
    benchmark:
        f"{OUT_DIR}/benchmark/01_download_genome.txt"
    shell:
        """
        mkdir -p {DATA_DIR}

        echo "Downloading genome from {params.genome_url}..." >> {log}
        wget -q -O {output.genome} {params.genome_url} 2>> {log} || \
        (echo "Error downloading genome" >> {log} && exit 1)
        
        echo "Downloading assembly report from {params.assembly_report_url}..." >> {log}
        wget -q -O {output.assembly_report} {params.assembly_report_url} 2>> {log} || \
        (echo "Error downloading assembly report" >> {log} && exit 1)
        
        echo "Download complete." >> {log}
        """

rule unzip_genome:
    """
    Unzip downloaded genome.
    """
    input:
        genome = f"{DATA_DIR}/{GENOME_FILE}"
    output:
        genome = f"{DATA_DIR}/{GENOME_FILE[:-3]}"
    log:
        f"{OUT_DIR}/log/02_unzip_genome.log"
    benchmark:
        f"{OUT_DIR}/benchmark/02_unzip_genome.txt"
    shell:
        """
        echo "Unzipping downloaded genome..." >> {log}
        gunzip {input.genome} 2>> {log} || \
        (echo "Error unzipping downloaded genome" >> {log} && exit 1)

        echo "Unzip complete." >> {log}
        """

rule subset_genome:
    """
    Select primary assembly and mitochondrial chromosomes.
    """
    input:
        genome = f"{DATA_DIR}/{GENOME_FILE[:-3]}",
        assembly_report = f"{DATA_DIR}/{ASSEMBLY_FILE}"
    output:
        renamed = f"{DATA_DIR}/genome_alignment.fa",
        ids = temp(f"{DATA_DIR}/subset_ids.txt"),
        genome = temp(f"{DATA_DIR}/genome_subset.fa"),
    log:
        f"{OUT_DIR}/log/03_subset_genome.log"
    benchmark:
        f"{OUT_DIR}/benchmark/03_subset_genome.txt"
    container:
        "docker://nottuh/sed-samtools:1.21"
    shell:
        """
        echo "Step 1: Extracting subset sequence IDs from assembly report..." >> {log}
        sort -k1,1V {input.assembly_report} |
        awk -F "\\t" '$8 == "Primary Assembly" || $8 == "non-nuclear" {{print $7}}' > {output.ids} 2>> {log} || \
        (echo "Error extracting IDs" >> {log} && exit 1)

        echo "Step 2: Extracting genome subset using samtools..." >> {log}
        samtools faidx {input.genome} -r {output.ids} -o {output.genome} 2>> {log} || \
        (echo "Error during genome subset extraction" >> {log} && exit 1)

        echo "Step 3: Replacing FASTA headers with UCSC-style headers..." >> {log}
        awk -v FS="\\t" 'NR==FNR {{header[">"$7] = ">"$10; next}} $0 ~ "^>" {{sub($0, header[$0])}} 1' \
        {input.assembly_report} {output.genome} > {output.renamed} 2>> {log} || \
        (echo "Error replacing FASTA headers" >> {log} && exit 1)

        echo "Subset and renaming complete." >> {log}
        """

rule index_genome:
    """
    Index the subsetted genome.
    """
    input:
        genome = f"{DATA_DIR}/genome_alignment.fa"
    output:
        index = [f"{DATA_DIR}/genome_index.{suffix}" for suffix in index_suffixes]
    params:
        stem = f"{DATA_DIR}/genome_index"
    log:
        f"{OUT_DIR}/log/04_index_genome.log"
    benchmark:
        f"{OUT_DIR}/benchmark/04_index_genome.txt"
    container:
        "docker://nottuh/sed-bowtie2:2.5.4"
    shell:
        """
        echo "Index the genome fasta file..." >> {log}
        bowtie2-build {input.genome} {params.stem} 2>> {log} || \
        (echo "Error indexing genome" >> {log} && exit 1)

        echo "Indexing of genome complete." >> {log}
        """

# ------------------------------- #
# 02. ChIP-seq Data Processing    #
# ------------------------------- #

rule download_chipseq:
    """
    Download ChIP-Seq samples.
    """
    output:
        fastq = f"{DATA_DIR}/{{srr}}.fastq"
    params:
        srr = lambda wildcards: wildcards.srr
    log:
        f"{OUT_DIR}/log/05_download_chipseq_{{srr}}.log"
    benchmark:
        f"{OUT_DIR}/benchmark/05_download_chipseq_{{srr}}.txt"
    container:
        "docker://nottuh/sed-sratools:3.2.1"
    shell:
        """
        echo "Downloading SRR file: {params.srr}..." >> {log}
        prefetch {params.srr} 2>> {log} || \
        (echo "Error running prefetch" >> {log} && exit 1)
        fastq-dump {params.srr} -O {DATA_DIR} 2>> {log} || \
        (echo "Error running fastq-dump" >> {log} && exit 1)

        echo "Download of {params.srr} complete." >> {log}
        """

rule align_reads:
    """
    Align reads to reference genome and processing them.
    """
    input:
        fastq = f"{DATA_DIR}/{{srr}}.fastq",
        index_files = [f"{DATA_DIR}/genome_index.{suffix}" for suffix in index_suffixes],
    output:
        sorted = temp(f"{DATA_DIR}/bams/{{srr}}_sorted.bam"),
        bam = f"{DATA_DIR}/bams/{{srr}}_markdup.bam",
        bai = f"{DATA_DIR}/bams/{{srr}}_markdup.bam.bai",
        sam = temp(f"{DATA_DIR}/bams/{{srr}}.sam"),
    params:
        index = f"{DATA_DIR}/genome_index",
    log:
        f"{OUT_DIR}/log/06_align_reads_{{srr}}.log"
    benchmark:
        f"{OUT_DIR}/benchmark/06_align_reads_{{srr}}.txt"
    container:
        "docker://nottuh/sed-bowtie2-samtools:latest"
    shell:
        """
        mkdir -p {DATA_DIR}/bams

        echo "Aligning ChIP-Seq reads for {wildcards.srr} to genome..." >> {log}
        bowtie2 -x {params.index} -U {input.fastq} -S {output.sam} 2>> {log} || \
        (echo "Error in alignment" >> {log} && exit 1)

        echo "Converting .sam to sorted .bam..." >> {log}
        samtools sort {output.sam} -o {output.sorted} 2>> {log} || \
        (echo "Error sorting .bam" >> {log} && exit 1)

        echo "Removing PCR duplicates..." >> {log}
        samtools markdup -r {output.sorted} {output.bam} 2>> {log} || \
        (echo "Error removing duplicates" >> {log} && exit 1)

        echo "Indexing .bam file..." >> {log}
        samtools index {output.bam} 2>> {log} || \
        (echo "Error indexing .bam" >> {log} && exit 1)

        echo "Processing of {wildcards.srr} complete." >> {log}
        """

rule move_bams:
    """
    Move non-control .bam and .bam.bai files.
    """
    input:
        expand(f"{DATA_DIR}/bams/{{srr}}_markdup.bam", srr=TARGET_IDS) +
        expand(f"{DATA_DIR}/bams/{{srr}}_markdup.bam.bai", srr=TARGET_IDS),
    output:
        expand(f"{DATA_DIR}/{{srr}}_markdup.bam", srr=TARGET_IDS) +
        expand(f"{DATA_DIR}/{{srr}}_markdup.bam.bai", srr=TARGET_IDS),
    log:
        f"{OUT_DIR}/log/07_move_bams.log"
    benchmark:
        f"{OUT_DIR}/benchmark/07_move_bams.txt"
    run:
        from pathlib import Path

        data_dir = Path(DATA_DIR)
        log_file = Path(log[0])

        with open(log_file, "w") as f:
            for srr in TARGET_IDS:
                srr_bam = data_dir / "bams" / f"{srr}_markdup.bam"
                srr_bai = data_dir / "bams" / f"{srr}_markdup.bam.bai"

                new_bam = data_dir / f"{srr}_markdup.bam"
                new_bai = data_dir / f"{srr}_markdup.bam.bai"

                f.write(f"Moving {srr_bam} to {new_bam}...")
                srr_bam.rename(new_bam)

                f.write(f"Moving {srr_bai} to {new_bai}...")
                srr_bai.rename(new_bai)

            f.write("All files moved successfully.\n")

rule call_peaks:
    """
    Call peaks with MACS2 using the target BAM files.
    """
    input:
        bam = f"{DATA_DIR}/{{srr}}_markdup.bam",
        bai = f"{DATA_DIR}/{{srr}}_markdup.bam.bai",
        control = lambda wildcards: f"{DATA_DIR}/bams/{CONTROLS_IDS[0]}_markdup.bam"
    output:
        narrowpeak = f"{DATA_DIR}/macs2/{{srr}}_peaks.broadPeak"
    params:
        broad_cutoff = BROAD_CUTOFF,
        extension_size = EXTENSION_SIZE,
        genome_size = GENOME_SIZE,
        prefix = lambda wildcards: wildcards.srr,
    log:
        f"{OUT_DIR}/log/08_call_peaks_{{srr}}.log"
    benchmark:
        f"{OUT_DIR}/benchmark/08_call_peaks_{{srr}}.txt"
    container:
        "docker://nottuh/macs2:2.2.7.1"     # change to nottuh/sed-macs2:2.2.7.1
    shell:
        """
        mkdir -p {DATA_DIR}/macs2

        echo "Calling peaks for {wildcards.srr}..." >> {log}
        macs2 callpeak \
            -t {input.bam} \
            -c {input.control} \
            -f BAM \
            -g {params.genome_size} \
            -n {params.prefix} \
            --broad \
            --broad-cutoff {params.broad_cutoff} \
            --outdir {DATA_DIR}/macs2 \
            --nomodel \
            --extsize {params.extension_size} \
            --keep-dup all 2>> {log} || \
        (echo "MACS2 peak calling failed" >> {log} && exit 1)

        echo "Peak calling for {wildcards.srr} complete." >> {log}
        """

rule convert_broadpeak:
    """
    Convert MACS2 broadPeak file to GFF3 format.
    """
    input:
        broadpeak = f"{DATA_DIR}/macs2/{{srr}}_peaks.broadPeak"
    output:
        gff3 = f"{DATA_DIR}/{{srr}}_constituent_enhancers.gff3"
    log:
        f"{OUT_DIR}/log/convert_broadpeak_{{srr}}.log"
    benchmark:
        f"{OUT_DIR}/benchmark/09_convert_broadpeak_{{srr}}.txt"
    shell:
        """
        echo "Converting broadPeak file for {wildcards.srr}..." >> {log}
        awk 'BEGIN{{OFS="\\t"}} {{print $1, $4, ".", $2+1, $3, ".", ".", ".", $4}}' {input.broadpeak} > {output.gff3} 2>> {log} || \
        (echo "broadPeak file conversion failed" >> {log} && exit 1)

        echo "Conversion of ROSE-compatible GFF3 file complete." >> {log}
        """
