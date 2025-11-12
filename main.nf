#!/usr/bin/env nextflow

include { call_all ; call_snps ; minimap2 ; get_clair3_model ; clair3 } from './process/variant_calling.nf'
include { apply_filters } from './process/consensus.nf'
include { apply_filters as apply_filters_clair3 } from './process/consensus.nf'
include { make_consensus ; make_clair3_consensus } from './process/consensus.nf'


workflow {
    if (params.help) {
        log.info(
            """
            ========================================================================
            Sundial Workflow

            Parameters:
            ------------------------------------------------------------------------
            --workflow  Workflow to run: bcftools or clair3
            --input_dir    Path to the sample's directory
            --ref_fasta   Path to the reference fasta
            --clair3_models_dir   Path to folder containing clair3 models
            --basecalling_model   model used for basecalling if clair3 should be used

            """.stripIndent()
        )
        exit(0)
    }

    input_files = Channel
        .fromPath("${params.input_dir}/${params.input_single_suffix}", checkIfExists: true)
        .ifEmpty { error("cannot find any reads matching ${params.input_single_suffix} in ${params.input_dir}") }
        .map { it -> tuple(it.getName().replaceFirst(/(?i)\.(fastq|fq)\.gz$/, ""), it) }

    input_files.take(3).view()

    read_ch = input_files.map { it ->
        [
            it[0],
            it[1],
            params.reference,
            params.accession,
            file(params.ref_fasta),
        ]
    }

    clair3_models_dir = params.clair3_models_dir.startsWith("/")
        ? params.clair3_models_dir
        : "${projectDir}/${params.clair3_models_dir}"

    if (params.workflow == "clair3") {
        rundial(read_ch, clair3_models_dir, params.basecalling_model)
    }
    else if (params.workflow == "bcftools") {
        rundial_with_bcftools(read_ch)
    }
    else {
        error("Please provide a valid workflow: bcftools or clair3")
    }
}

// Main workflow uses clair3 variant calling with bcftools to fill in the gaps
workflow rundial {
    take:
    fastq_files
    clair3_models_dir
    dorado_model

    main:

    ref_name_ch = fastq_files.map { it -> [it[0], it[2]] }

    fastq_files = fastq_files.map { it -> tuple(it[0], it[1], it[4]) }

    clair3_model = get_clair3_model(dorado_model)

    valid_model_ch = fastq_files
        .combine(clair3_model)
        .filter { it -> it[3] != "failed" }
        .map { it -> tuple(it[0], it[1], it[2]) }
    invalid_model_ch = fastq_files
        .combine(clair3_model)
        .filter { it -> it[3] == "failed" }
        .map { it -> tuple(it[0], it[1], it[2]) }

    valid_model_ch.take(1).view { "Appropriate Clair3 Model Found" }
    invalid_model_ch.take(1).view { "No appropriate Clair3 Model Found" }

    // If no appropriate clair3 model is found, use bcftools
    rundial_with_bcftools(invalid_model_ch)

    // If appropriate clair3 model is found, run clair3
    model_path = clair3_model.map { it -> clair3_models_dir + "/${it}.tar.gz" }
    model_path.view { "Using Clair3 Model: ${it}" }

    minimap2(valid_model_ch)
    clair3(minimap2.out.sorted_alignment, model_path)

    clair3_filter_params = Channel.fromPath("${moduleDir}/process/clair3_filter_params.yml").first()
    apply_filters_clair3(clair3.out.vcf, clair3_filter_params, "clair3_")

    // Still use normal params with bcftools support vcf
    filter_params = Channel.fromPath("${moduleDir}/process/filter_params.yml").first()
    call_snps(minimap2.out.sorted_alignment)
    apply_filters(call_snps.out.gvcf, filter_params, "bcftools_")


    consensus_params = Channel.fromPath("${moduleDir}/process/clair3_consensus_params.yml").first()
    calls = apply_filters.out.filtered_gvcf.join(apply_filters_clair3.out.filtered_gvcf)
    make_clair3_consensus(calls, consensus_params)

    emit:
    alignment = rundial_with_bcftools.out.alignment.concat(minimap2.out.sorted_alignment).join(ref_name_ch).map { it -> tuple(it[0], it[1], it[3], null, it[2]) }
    gvcf = rundial_with_bcftools.out.gvcf.concat(clair3.out.vcf).join(ref_name_ch)
    final_fasta = rundial_with_bcftools.out.final_fasta.concat(make_clair3_consensus.out.final_fasta).join(ref_name_ch)
    variable_length_fasta = rundial_with_bcftools.out.variable_length_fasta.concat(make_clair3_consensus.out.variable_length_fasta).join(ref_name_ch)
    final_vcf = rundial_with_bcftools.out.final_vcf.concat(make_clair3_consensus.out.final_vcf).join(ref_name_ch)
    full_vcf = rundial_with_bcftools.out.full_vcf.concat(make_clair3_consensus.out.full_vcf).join(ref_name_ch)
    creation_report_json = rundial_with_bcftools.out.creation_report_json.concat(make_clair3_consensus.out.report_json).join(ref_name_ch)
}


workflow rundial_with_bcftools {
    take:
    fastq_files

    main:
    filter_params = Channel.fromPath("${moduleDir}/process/filter_params.yml").first()
    consensus_params = Channel.fromPath("${moduleDir}/process/consensus_params.yml").first()

    minimap2(fastq_files)
    call_all(minimap2.out.sorted_alignment)
    apply_filters(call_all.out.gvcf, filter_params, "")

    calls = apply_filters.out.filtered_gvcf
    calls = calls.join(fastq_files).map { it -> tuple(it[0], it[1], it[3]) }
    calls.view { "Calls channel: ${it}" }
    make_consensus(calls, consensus_params)

    emit:
    alignment = minimap2.out.sorted_alignment
    gvcf = apply_filters.out.filtered_gvcf
    final_fasta = make_consensus.out.final_fasta
    variable_length_fasta = make_consensus.out.variable_length_fasta
    final_vcf = make_consensus.out.final_vcf
    full_vcf = make_consensus.out.full_vcf
    creation_report_json = make_consensus.out.report_json
}
