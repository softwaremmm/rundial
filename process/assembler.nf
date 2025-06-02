params.test_cpus = ""

process minimap2 {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    cpus {
        params.testing == "" ? 10 : (params.test_cpus == "" ? 2 : params.test_cpus)
    }
    container {
        params.test_container_rundial == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/rundial:221bf51' : params.test_container_rundial
    }

    pod label: "name", value: "rundial:minimap2"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(fq)
    path reference

    output:
    tuple val(sample_name), path("final.bam"), emit: sorted_alignment

    script:
    """
    minimap2 -t ${task.cpus} -a -L --sam-hit-only --secondary=no -x map-ont ${reference} ${fq} | \
        samtools sort -@ ${task.cpus} -o final.bam
    """
}

process call_snps {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    cpus {
        params.testing == "" ? 10 : (params.test_cpus == "" ? 2 : params.test_cpus)
    }
    container {
        params.test_container_rundial == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/rundial:221bf51' : params.test_container_rundial
    }

    pod label: "name", value: "rundial:call_snps"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(sorted_alignment)
    path reference

    output:
    tuple val(sample_name), path("calls.gvcf.gz"), emit: gvcf
    tuple val(sample_name), path("calls.raw_variants.vcf.gz"), emit: variants

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

    bcftools view -v snps calls.gvcf.gz -Oz -o calls.raw_variants.vcf.gz

    # Clean up large files to save space locally
    rm -r pileups
    rm pileup.bcf
    """
}

process call_all {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "." + filename }
    cpus {
        params.testing == "" ? 10 : (params.test_cpus == "" ? 2 : params.test_cpus)
    }
    container {
        params.test_container_rundial == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/rundial:221bf51' : params.test_container_rundial
    }

    pod label: "name", value: "rundial:call_all"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(sorted_alignment)
    path reference

    output:
    tuple val(sample_name), path("calls.gvcf.gz"), emit: gvcf
    tuple val(sample_name), path("calls.raw_variants.vcf.gz"), emit: variants

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

    bcftools view -v snps,indels calls.gvcf -o calls.raw_variants.vcf

    # Clean up large files to save space locally
    rm -r pileups
    rm pileup.bcf
    rm snps.gvcf.gz
    rm indels.gvcf.gz
    bgzip calls.gvcf
    bgzip calls.raw_variants.vcf
    """
}


process get_clair3_model {
    container {
        params.test_container_rundial == "" ? 'lhr.ocir.io/lrbvkel2wjot/gpas/rundial:221bf51' : params.test_container_rundial
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
    cpus {
        params.testing == "" ? 4 : (params.test_cpus == "" ? 2 : params.test_cpus)
    }
    container "hkubal/clair3:v1.0.5"

    pod label: "name", value: "rundial:clair3"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path(bam)
    path reference
    path model

    output:
    tuple val(sample_name), path("clair3.vcf.gz"), emit: vcf

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
        --output=clair_out > clair3.log

    mv clair_out/merge_output.vcf.gz clair3.vcf.gz
    """
}
