configfile: "config/example.yaml"


from pathlib import Path


# Initialising data variables
assembly = config["data"]["assembly"]
genome = config["data"]["genome"]

rule all:
    input:
        f"data/{Path(assembly).name}",
        f"data/{Path(genome).stem}",
        "data/genome_alignment.fa",
 
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
        awk -v FS="\t" 'NR==FNR {{header[">"$7] = ">"$10; next}} $0 ~ "^>" {{sub($0, header[$0]); print}}1' {input.assembly} {params.genome} > {output} 2>&1 | tee -a {log}
        """
