process download_data {

    publishDir path: "${params.outdir}/data", mode: "copy"

    input:
        val(file)

    output:
        path("*")

    script:
        ext = file(file).extension
        """
        wget -q ${file}

        if [[ ${ext} == "gz" ]]; then
            gunzip *.gz
        fi
        """
}