process BWA_MEM2_INDEX {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bwa-mem2=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa-mem2:2.2.1--he513fc3_0':
        'biocontainers/bwa-mem2:2.2.1--he513fc3_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("${fasta}.*"), emit: index
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bwa-mem2 index ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version | head -n1 | sed 's/bwa-mem2 //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.0123
    touch ${fasta}.amb
    touch ${fasta}.ann
    touch ${fasta}.bwt.2bit.64
    touch ${fasta}.pac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version | head -n1 | sed 's/bwa-mem2 //g')
    END_VERSIONS
    """
}
