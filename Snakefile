from pathlib import Path


configfile: "config/example.yaml"


# Get data from configuration file
DATA_DIR = config["dir"]["data_dir"]
OUT_DIR = config["dir"]["output_dir"]

# Initialising variables
GENOME_FILE = Path(config["data"]["genome"]).name
ASSEMBLY_FILE = Path(config["data"]["assembly_report"]).name


rule all:
    """
    Defining the final expected output files.
    """
    input:
        f"{DATA_DIR}/{GENOME_FILE[:-3]}",
        f"{DATA_DIR}/{ASSEMBLY_FILE}", 
        f"{DATA_DIR}/GRCh38_alignment.fa",


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
        f"{OUT_DIR}/log/download_genome.log"
    benchmark:
        f"{OUT_DIR}/benchmark/download_genome.txt"
    shell:
        """
        mkdir -p {DATA_DIR}

        echo "Downloading genome from {params.genome_url}" >> {log}
        wget -q -O {output.genome} {params.genome_url} 2>> {log} || (echo "Error downloading genome" >> {log} && exit 1)
        
        echo "Downloading assembly report from {params.assembly_report_url}" >> {log}
        wget -q -O {output.assembly_report} {params.assembly_report_url} 2>> {log} || (echo "Error downloading assembly report" >> {log} && exit 1)
        
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
        f"{OUT_DIR}/log/unzip_genome.log"
    benchmark:
        f"{OUT_DIR}/benchmark/unzip_genome.txt"
    shell:
        """
        echo "Unzipping downloaded genome" >> {log}
        gunzip {input.genome} 2>> {log} || (echo "Error unzipping downloaded genome" >> {log} && exit 1)

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
        renamed = f"{DATA_DIR}/GRCh38_alignment.fa"
    params:
        ids = temp(f"{DATA_DIR}/subset_ids.txt"),
        genome = temp(f"{DATA_DIR}/genome_subset.fa")
    log:
        f"{OUT_DIR}/log/subset_genome.log"
    benchmark:
        f"{OUT_DIR}/benchmark/subset_genome.txt"
    container:
        "docker://nottuh/sed-samtools:1.21"
    shell:
        """
        echo "Step 1: Extracting subset sequence IDs from assembly report..." >> {log}
        sort -k1,1V {input.assembly_report} |
        awk -F "\\t" '$8 == "Primary Assembly" || $8 == "non-nuclear" {{print $7}}' > {params.ids} 2>> {log} || \
        (echo "Error extracting IDs" >> {log} && exit 1)

        echo "Step 2: Extracting genome subset using samtools..." >> {log}
        samtools faidx {input.genome} -r {params.ids} -o {params.genome} 2>> {log} || \
        (echo "Error during genome subset extraction" >> {log} && exit 1)

        echo "Step 3: Replacing FASTA headers with UCSC-style headers..." >> {log}
        awk -v FS="\\t" 'NR==FNR {{header[">"$7] = ">"$10; next}} $0 ~ "^>" {{sub($0, header[$0]); print}}1' \
        {input.assembly_report} {params.genome} > {output.renamed} 2>> {log} || \
        (echo "Error replacing FASTA headers" >> {log} && exit 1)

        echo "Subset and renaming complete." >> {log}
        """