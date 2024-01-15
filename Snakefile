configfile: "config/example.yaml"


from pathlib import Path


# Initialising data variables
assembly = config["data"]["assembly"]
genome = config["data"]["genome"]

rule all:
    input:
        f"data/{Path(assembly).name}",
        f"data/{Path(genome).stem}"
 
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
