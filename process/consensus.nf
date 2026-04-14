process apply_filters {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? params.container_prefix + '/gpas/rundial:0c50b02' : params.test_container_rundial
    }
    cpus 1

    pod label: "name", value: "rundial:apply_filters"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(gvcf)
    path filter_params
    val caller

    output:
    tuple val(sample_name), path("alternate-${caller}.vcf.gz"), emit: filtered_vcf

    script:
    """
    rundial filter --verbose --overwrite \
        -p ${filter_params} \
        -o alternate-${caller}.vcf -i ${gvcf}
    bgzip alternate-${caller}.vcf
    """
}

process make_consensus {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? params.container_prefix + '/gpas/rundial:0c50b02' : params.test_container_rundial
    }
    cpus 1
    memory "4 GB"

    pod label: "name", value: "rundial:make_consensus"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(alternate_bcftools_vcf)
    path reference
    path consensus_params

    output:
    tuple val(sample_name), path("final.fasta"), emit: final_fasta
    tuple val(sample_name), path("final.variable_length.fasta"), emit: variable_length_fasta
    tuple val(sample_name), path("variants.vcf"), emit: variants_vcf
    tuple val(sample_name), path("all_calls.vcf"), emit: all_calls_vcf
    tuple val(sample_name), path("genome_creation_report.json"), emit: report_json

    script:
    """
    rundial consensus --verbose \
        -p ${consensus_params} \
        --ref-fasta ${reference} \
        -i ${alternate_bcftools_vcf} \
        -o consensus

    bcftools view -v snps,indels consensus.vcf > variants.vcf

    # Rename files to match output specification
    mv consensus.vcf all_calls.vcf
    mv consensus.report.json genome_creation_report.json
    mv consensus.fasta final.fasta
    mv consensus.variable_length.fasta final.variable_length.fasta

    # replace header of fasta files
    sed -i "s/^>/>${sample_name}:/" final.fasta
    sed -i "s/^>/>${sample_name}:/" final.variable_length.fasta
    """
}

process make_clair3_consensus {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? params.container_prefix + '/gpas/rundial:0c50b02' : params.test_container_rundial
    }
    cpus 1
    memory "3 GB"

    pod label: "name", value: "rundial:make_clair3_consensus"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path("alternate-bcftools.vcf.gz"), path("clair3.vcf.gz")
    path reference
    path consensus_params

    output:
    tuple val(sample_name), path("final.fasta"), emit: final_fasta
    tuple val(sample_name), path("final.variable_length.fasta"), emit: variable_length_fasta
    tuple val(sample_name), path("variants.vcf"), emit: variants_vcf
    tuple val(sample_name), path("all_calls.vcf"), emit: all_calls_vcf
    tuple val(sample_name), path("genome_creation_report.json"), emit: report_json

    script:
    """
    rundial consensus --verbose \
        -p ${consensus_params} \
        --ref-fasta ${reference} \
        -i clair3.vcf.gz \
        -s alternate-bcftools.vcf.gz \
        -o consensus

    bcftools view -v snps,indels consensus.vcf > variants.vcf

    # Rename files to match output specification
    mv consensus.vcf all_calls.vcf
    mv consensus.report.json genome_creation_report.json
    mv consensus.fasta final.fasta
    mv consensus.variable_length.fasta final.variable_length.fasta

    # replace header of fasta files
    sed -i "s/^>/>${sample_name}:/" final.fasta
    sed -i "s/^>/>${sample_name}:/" final.variable_length.fasta
    """
}


process reassess_genome_creation {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? params.container_prefix + '/gpas/rundial:0c50b02' : params.test_container_rundial
    }
    cpus 1

    pod label: "name", value: "rundial:reassess_genome_creation"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fasta), path(vcf), path(mask)
    val min_allele_dp
    val min_allele_pc
    val min_strand_pc
    val cluster_window

    output:
    tuple val(sample_name), path("updated_creation_report.json"), emit: report_json

    script:
    """
    assess_genome_creation --fasta ${fasta} --vcf ${vcf} \
        --mask ${mask} \
        --min_allele_dp ${min_allele_dp} \
        --min_allele_pc ${min_allele_pc} \
        --min_strand_pc ${min_strand_pc} \
        --cluster_window ${cluster_window} \
        -o updated_creation_report.json
    """
}
