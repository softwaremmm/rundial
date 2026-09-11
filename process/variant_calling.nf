params.test_cpus = ""

process minimap2 {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? params.container_prefix + '/gpas/rundial:0f42216' : params.test_container_rundial
    }
    cpus {
        params.testing == "" ? 4 : params.test_cpus
    }

    pod label: "name", value: "rundial:minimap2"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), val(ref_id), path(reference), path(fq)

    output:
    tuple val(sample_name), val(ref_id), path(reference), path("final.bam"), emit: sorted_alignment

    script:
    """
    set -o pipefail

    minimap2 -t ${task.cpus} -a -L --sam-hit-only --secondary=no -x map-ont ${reference} ${fq} | \
        samtools sort -@ ${task.cpus} -o final.bam
    """
}

process call_snps {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? params.container_prefix + '/gpas/rundial:0f42216' : params.test_container_rundial
    }
    cpus {
        params.testing == "" ? 4 : params.test_cpus
    }

    pod label: "name", value: "rundial:call_snps"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), val(ref_id), path(reference), path(sorted_alignment)

    output:
    tuple val(sample_name), val(ref_id), path(reference), path("calls.gvcf.gz"), emit: gvcf

    script:
    """
    echo "Running bcftools mpileup in parallel"
    date +"%T"

    # script also makes indexes
    multithreaded_bcftools -a ${sorted_alignment} -r ${reference} \
        -o pileup.bcf -t ${task.cpus} \
        --settings "-I -x -Q 10 -a INFO/SCR,INFO/ADR,INFO/ADF,FORMAT/SP,FORMAT/AD -h100 -M10000"

    echo "Running bcftools call"
    date +"%T"
    # Report all alleles during calling
    bcftools call --threads ${task.cpus} --ploidy 1 -m -A -V indels pileup.bcf \
        | bcftools filter -e 'DP==0' \
        | bcftools norm -d exact -f ${reference} -Oz -o calls.gvcf.gz

    echo "Finished Calling Variants"
    date +"%T"

    # Clean up large files to save space locally
    rm -r pileups
    rm pileup.bcf
    """
}

process call_all {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? params.container_prefix + '/gpas/rundial:0f42216' : params.test_container_rundial
    }
    cpus {
        params.testing == "" ? 4 : params.test_cpus
    }

    pod label: "name", value: "rundial:call_all"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), val(ref_id), path(reference), path(sorted_alignment)

    output:
    tuple val(sample_name), val(ref_id), path(reference), path("calls.gvcf.gz"), emit: gvcf

    script:
    """
    echo "Running bcftools mpileup in parallel"

    # script also makes indexes
    multithreaded_bcftools -a ${sorted_alignment} -r ${reference} \
        -o pileup.bcf -t ${task.cpus} \
        --settings "-x -Q 10 -a INFO/SCR,INFO/ADR,INFO/ADF,FORMAT/SP,FORMAT/AD -h100 -M10000"

    echo "Running bcftools call"
    date +"%T"
    # Report all alleles during calling
    bcftools call --threads ${task.cpus} --ploidy 1 -m -A -V indels pileup.bcf \
        | bcftools filter -e 'DP==0' \
        | bcftools norm -d exact -f ${reference} -Oz -o snps.gvcf.gz
    bcftools call --threads ${task.cpus} --ploidy 1 -m -A -v -V snps pileup.bcf \
        | bcftools filter -e 'DP==0' \
        | bcftools norm -d exact -f ${reference} -Oz -o indels.gvcf.gz
    bcftools index snps.gvcf.gz
    bcftools index indels.gvcf.gz
    bcftools concat -a snps.gvcf.gz indels.gvcf.gz -o calls.gvcf

    echo "Finished Calling Variants"
    date +"%T"

    # Clean up large files to save space locally
    rm -r pileups
    rm pileup.bcf
    rm snps.gvcf.gz
    rm indels.gvcf.gz
    bgzip calls.gvcf
    """
}


process get_clair3_model {
    container {
        params.test_container_rundial == "" ? params.container_prefix + '/gpas/rundial:0f42216' : params.test_container_rundial
    }

    pod label: "name", value: "rundial:get_clair3_model"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    val dorado_model

    output:
    stdout

    script:
    """
    clair3_model=\$(rundial dorado_to_clair3_model ${dorado_model})
    echo -n \$clair3_model
    """
}

process clair3 {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    container {
        params.test_container_rundial == "" ? params.container_prefix + '/gpas/rundial:0f42216' : params.test_container_rundial
    }
    cpus {
        params.testing == "" ? 4 : params.test_cpus
    }

    pod label: "name", value: "rundial:clair3"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), val(ref_id), path(reference), path(bam)
    path model

    output:
    tuple val(sample_name), val(ref_id), path(reference), path("clair3.vcf.gz"), emit: vcf

    script:
    """
    samtools index ${bam}
    samtools faidx ${reference}

    tar -xzf ${model}
    model_folder=\$(basename ${model} .tar.gz)

    run_clair3.sh \
        --bam_fn=${bam} \
        --ref_fn=${reference} \
        --threads=${task.cpus} \
        --platform="ont" \
        --model_path=\$model_folder \
        --include_all_ctgs \
        --no_phasing_for_fa \
        --haploid_sensitive \
        --enable_long_indel \
        --qual=2 \
        --print_ref_calls \
        --sample_name=${bam} \
        --output=clair_out

    mv clair_out/merge_output.vcf.gz clair3.vcf.gz
    """
}
