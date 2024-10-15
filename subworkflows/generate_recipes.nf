#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process subset_reference_accessions {
    label "process_low"

    container 'biocontainers/python:3.10.4'

    input:
    path reference_csv
    val sample_size

    output:
    path "references.csv"

    script:
    """
    subset_accessions.py ${reference_csv} "references.csv" category_id ${sample_size}
    """
}

process subset_dataset_accessions {
    label "process_low"

    container 'biocontainers/python:3.10.4'

    input:
    path dataset_csv
    val sample_size

    output:
    path "dataset.csv"

    script:
    """
    subset_accessions.py ${dataset_csv} "dataset.csv" platform ${sample_size}
    """
}

process download_reference_fasta {
    label "process_low"

    container "community.wave.seqera.io/library/ncbi-datasets-cli_unzip:ec913708564558ae"

    input:
    tuple val(accession), val(category)    

    output:
    tuple val(accession), val(category), path("${accession}_genomic.fna")

    
    script:

    """
    datasets download genome accession ${accession}
    unzip -o ncbi_dataset.zip
    output="\$(readlink -f ncbi_dataset/data/*/*_genomic.fna)"
    mv "\${output}" ${accession}_genomic.fna
    """

}

process download_dataset_accession {

    container "biocontainers/sra-tools:2.7.0--0"

    storeDir "${params.dataset_dir}/"
    input:
    tuple val(accession), val(platform)

    output:
    tuple val(accession), val(platform), path("${accession}/*_pass.fastq.gz")

    script:
    """
    prefetch ${accession}
    fastq-dump --outdir ${accession} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip */*.sra
    """
}

process download_dataset_accession_paired {
    label "process_low"
    
    container "biocontainers/sra-tools:2.7.0--0"

    storeDir "${params.dataset_dir}/"
    input:
    tuple val(accession), val(platform)

    output:
    tuple val(accession), val(platform), path("${accession}/*_pass_1.fastq.gz"), path("${accession}/*_pass_2.fastq.gz")

    script:
    """
    prefetch ${accession}
    fastq-dump --outdir ${accession} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip */*.sra
    """
}

workflow get_base_fastq {
    take:
        dataset_row
    main:
        dataset_row.branch { accession, platform, index, reads1 ->
            download: reads1 == ""
                return tuple(accession, platform, index)
            other: true
                return tuple(index, accession, platform, reads1)
            }.set { result }
        result.download.tap{ to_download }
        download_dataset_accession(to_download.map{ accession, platform, index -> [accession, platform] }.unique())
        result.download.combine(download_dataset_accession.out, by: 0).map{ accession, platform, index, platform1, reads1 -> [index, accession, platform, reads1]}.set{downloaded}
        result.other.concat(downloaded).set{ base_fastq }
    emit:
        base_fastq
}

workflow get_base_fastq_paired {
    take:
        dataset_row
    main:
        dataset_row.branch { accession, platform, index, reads1, reads2 ->
            download: reads1 == ""
                return tuple(accession, platform, index)
            other: true
                return tuple(index, accession, platform, reads1, reads2)
            }.set { result }
        result.download.tap{ to_download }
        download_dataset_accession_paired(to_download.map{ accession, platform, index -> [accession, platform] }.unique())
        result.download.combine(download_dataset_accession_paired.out, by: 0).map{ accession, platform, index, platform1, reads1, reads2 -> [index, accession, platform, reads1, reads2]}.set{downloaded}
        result.other.concat(downloaded).set{ base_fastq }
    emit:
        base_fastq
}

workflow get_reference_fastas {
    take:
        reference_csv
    main:
        subset_reference_accessions(reference_csv, params.num_iterations)
        subset_reference_accessions.out.splitCsv(header: true).map { 
            row -> def file_path = row.path ? file(row.path) : null
                tuple("${row.accession}", "${row.category_id}", "${row.index}", file_path)
            }.set{ reference_accessions }

        reference_accessions.branch {
            // if file_path is not null, then that means the user supplied a path to the fasta, so branch into local
            local: it[3] != null
            // if it is null, then try and download based on accession
            download: it[3] == null
        }.set { reference }

        download_reference_fasta(reference.download.map{ accession, category, index, file_path -> [accession, category] }.unique())
        
        reference.download
            .combine(download_reference_fasta.out, by: 0)
            .map{ accession, category, index, file_path, category1, fasta -> [index, accession, category, fasta]}
            .set{ downloaded }

        reference.local
            .map{ accession, category, index, file_path -> [index, accession, category, file_path]}
            .mix(downloaded)
            .set{ all_references }
    emit:
      all_references
}

workflow get_base_datasets {
    take:
        dataset_csv
    main:
        subset_dataset_accessions(dataset_csv, params.num_iterations)
        subset_dataset_accessions.out.splitCsv(header: true).map{row -> ["${row.public_database_accession}","${row.platform}","${row.index}","${row.human_filtered_reads_1}","${row.human_filtered_reads_2}"]}.set{ dataset_accessions }

        dataset_accessions.branch { accession, platform, index, reads1, reads2 ->
            paired: platform == "illumina"
                return tuple(accession, platform, index, reads1, reads2)
            unpaired: true
                return tuple(accession, platform, index, reads1)
            }.set { by_platform }

        get_base_fastq_paired(by_platform.paired)
        get_base_fastq(by_platform.unpaired)

    emit:
        paired = get_base_fastq_paired.out
        unpaired = get_base_fastq.out
}

workflow generate_recipes {
    take:
        reference_csv
        dataset_csv
    main:
        get_reference_fastas(reference_csv)
        coverages = channel.from(params.coverages)
        get_reference_fastas.out.combine(coverages).set{ references }

        get_base_datasets(dataset_csv)
        references.combine(get_base_datasets.out.paired, by: 0).set{ paired_recipes }
        references.combine(get_base_datasets.out.unpaired, by: 0).set{ unpaired_recipes }
        paired_recipes.view()
        unpaired_recipes.view()
     emit:
        paired = paired_recipes
        unpaired = unpaired_recipes
}
