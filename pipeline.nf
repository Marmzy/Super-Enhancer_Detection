#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { download_data } from "./modules/preprocessing.nf"

workflow {

    // Download and unzip data
    ncbi_files = download_data(Channel.fromList([params.genome, params.assembly]))
    ncbi_files.view { "${it}" }

}