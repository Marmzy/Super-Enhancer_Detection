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


process index_genome {

    publishDir path: "${params.outdir}/data", mode: "copy"

    input:
        file(alignment)
        val(genome_basename)

    output:
        path("${genome_basename}_index*")

    script:
        """
        /home/calvin/bowtie2-2.2.2/bowtie2-build $alignment ${genome_basename}_index
        """
}


process replace_headers {

    publishDir path: "${params.outdir}/data", mode: "copy"

    input:
        file(assembly)
        file(subset_genome)
        val(genome_basename)

    output:
        path("*.fa")

    script:
        """
        awk -v FS="\t" 'NR==FNR {header[">"\$7] = ">"\$10; next} \$0 ~ "^>" {sub(\$0, header[\$0]); print}1' $assembly $subset_genome > ${genome_basename}_alignment.fa
        """
}


process subset_genome {

    publishDir path: "${params.outdir}/data", mode: "copy"

    input:
        tuple file(assembly), file(genome)

    output:
        path("subset_ids.txt")
        path("${genome_basename}_subset.fa"), emit: genome_subset
        val(genome_basename), emit: genome_basename

    script:
        genome_basename = genome.baseName
        """
        sort -k1,1V ${assembly} | awk -F "\t" '\$8 == "Primary Assembly" || \$8 == "non-nuclear" {print \$7}' > subset_ids.txt
        samtools faidx ${genome} -r subset_ids.txt -o ${genome_basename}_subset.fa
        """
}