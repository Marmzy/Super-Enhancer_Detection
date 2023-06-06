#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { download_data; subset_genome } from "./modules/preprocessing.nf"

workflow {

    // Downloading the Ch38 unpatched genome and assembly report
    ncbi_files = download_data(Channel.fromList([["genome", params.genome], ["assembly", params.assembly]]))
    ncbi_files = ncbi_files.toSortedList { a, b -> a[0] <=> b[0] } | map { it -> [it[0][1], it[1][1]]}

    // Select entries in the genome FASTA file, belonging to the Primary Assembly or Mitochondrial chromosome
    subsets = subset_genome(ncbi_files)

}