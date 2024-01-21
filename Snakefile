configfile: "config/example.yaml"


from pathlib import Path


# Initialising data variables
assembly = config["data"]["assembly"]
genome = config["data"]["genome"]

# Define wildcards
index_suffixes = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]

rule all:
    input:
        f"data/{Path(assembly).name}",
        f"data/{Path(genome).stem}",
        "data/genome_alignment.fa",
        expand("data/genome_index.{suffix}", suffix=index_suffixes),
 
rule download_assembly:
    output:
        "data/" + str(Path(assembly).name)
    log:
        "output/log/download_assembly.log"
    shell:
        """ 
        # Downloading the genome assembly report
        wget -O {output} {assembly} --tries=3 --waitretry=5 > {log} 2>&1 || echo "Error: Assembly download failed for {assembly}" >> {log}
        """

rule download_genome:
    output:
        "data/" + str(Path(genome).stem)
    log:
        "output/log/download_genome.log"
    shell:
        """
        # Downloading and unzipping the genome
        wget -O {output}.gz {genome} --tries=3 --waitretry=5 &&
        gunzip {output} 2>> {log}
        """

rule subset_genome:
    input:
        assembly = "data/" + str(Path(assembly).name),
        genome = "data/" + str(Path(genome).stem)
    output:
        "data/genome_alignment.fa"
    params:
        ids = temp("data/subset_ids.txt"),
        genome = temp("data/genome_subset.fa")
    log:
        "output/log/subset_genome.log"
    shell:
        """
        # Select the primary assembly accessions from the genome fasta file
        sort -k1,1V {input.assembly} | awk -F "\t" '$8 == "Primary Assembly" || $8 == "non-nuclear" {{print $7}}' > {params.ids} 2>&1 | tee -a {log}
        samtools faidx {input.genome} -r {params.ids} -o {params.genome} 2>&1 | tee -a {log}

        # Replacing the RefSeq-style fasta headers with UCSC-style headers
        awk -v FS="\t" 'NR==FNR {{header[">"$7] = ">"$10; next}} $0 ~ "^>" {{print header[$0]; next}}1' {input.assembly} {params.genome} > {output} 2>&1 | tee -a {log}
        """


rule index_genome:
    input:
        genome = "data/genome_alignment.fa"
    output:
        index = expand("data/genome_index.{suffix}", suffix=index_suffixes)
    params:
        stem = "data/genome_index"
    log:
        "output/log/index_genome.log"
    shell:
        """
        # Index the genome fasta file
        bowtie2-build {input.genome} {params.stem} 2>&1 | tee -a {log}
        """