process download_data {

    publishDir path: "${params.outdir}/data", mode: "copy"

    input:
        tuple val(name), val(file)

    output:
        tuple val(name), path("*")

    script:
        ext = file(file).extension
        """
        wget -q ${file}

        if [[ ${ext} == "gz" ]]; then
            gunzip *.gz
        fi
        """
}


process subset_genome {

    publishDir path: "${params.outdir}/data", mode: "copy"

    input:
        tuple file(assembly), file(genome)

    output:
        path("subset_ids.txt"), emit: subset_ids
        path("${genome_basename}_subset.fa"), emit: genome_subset

    script:
        genome_basename = genome.baseName
        """
        sort -k1,1V ${assembly} | awk -F "\t" '\$8 == "Primary Assembly" || \$8 == "non-nuclear" {print \$7}' > subset_ids.txt
        samtools faidx ${genome} -r subset_ids.txt -o ${genome_basename}_subset.fa
        """
}