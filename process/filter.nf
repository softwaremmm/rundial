process apply_filters {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + file_prefix + filename }
    container {
        params.test_container_rundial == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/rundial:8789e63' : params.test_container_rundial
    }

    pod label: "name", value: "rundial:apply_filters"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(gvcf)
    path filter_params
    val file_prefix

    output:
    tuple val(sample_name), path("filtered.gvcf.gz"), emit: filtered_gvcf

    script:
    """
    rundial filter --verbose --overwrite \
        -p ${filter_params} \
        -o filtered.gvcf -i ${gvcf}
    bgzip filtered.gvcf
    """
}

process make_consensus {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/rundial:8789e63' : params.test_container_rundial
    }

    pod label: "name", value: "rundial:make_consensus"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(filtered_gvcf)
    path reference
    path consensus_params

    output:
    tuple val(sample_name), path("final.full.fasta"), emit: full_consensus
    tuple val(sample_name), path("final.fasta"), emit: final_fasta
    tuple val(sample_name), path("final.indel.fasta"), emit: indel_fasta, optional: true
    tuple val(sample_name), path("final.vcf"), emit: final_vcf
    tuple val(sample_name), path("final.full.vcf"), emit: full_vcf
    tuple val(sample_name), path("genome_creation_report.json"), emit: report_json

    script:
    """
    rundial consensus --verbose \
        -p ${consensus_params} \
        --ref-fasta ${reference} \
        -i ${filtered_gvcf} \
        -o final

    mv final.vcf final.full.vcf
    bcftools view -v snps,indels -i 'GT!="0/0"' final.full.vcf > final.vcf
    mv final.report.json genome_creation_report.json
    mv final.variable_length.fasta final.indel.fasta

    # replace header of fasta files
    sed -i "s/^>/>${sample_name}:/" final.full.fasta
    sed -i "s/^>/>${sample_name}:/" final.fasta
    sed -i "s/^>/>${sample_name}:/" final.indel.fasta
    """
}

process make_clair3_consensus {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/rundial:8789e63' : params.test_container_rundial
    }

    pod label: "name", value: "rundial:make_clair3_consensus"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path("filtered.gvcf.gz"), path("clair3.vcf.gz")
    path reference
    path consensus_params

    output:
    tuple val(sample_name), path("final.full.fasta"), emit: full_consensus
    tuple val(sample_name), path("final.fasta"), emit: final_fasta
    tuple val(sample_name), path("final.indel.fasta"), emit: indel_fasta, optional: true
    tuple val(sample_name), path("final.vcf"), emit: final_vcf
    tuple val(sample_name), path("final.full.vcf"), emit: full_vcf
    tuple val(sample_name), path("genome_creation_report.json"), emit: report_json

    script:
    """
    rundial consensus --verbose \
        -p ${consensus_params} \
        --ref-fasta ${reference} \
        -i clair3.vcf.gz \
        -s filtered.gvcf.gz \
        -o final

    mv final.vcf final.full.vcf
    bcftools view -v snps,indels final.full.vcf > final.vcf
    mv final.report.json genome_creation_report.json
    mv final.variable_length.fasta final.indel.fasta

    # replace header of fasta files
    sed -i "s/^>/>${sample_name}:/" final.full.fasta
    sed -i "s/^>/>${sample_name}:/" final.fasta
    sed -i "s/^>/>${sample_name}:/" final.indel.fasta
    """
}
