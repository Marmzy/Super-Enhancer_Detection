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