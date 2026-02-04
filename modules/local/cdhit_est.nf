process CDHIT_EST {
    tag "Group $meta.id"
    label 'process_medium'

    conda "bioconda::cd-hit=4.8.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit%3A4.8.1--h5ca1c30_13':
        'biocontainers/cd-hit:4.8.1--h5ca1c30_13' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.dedup.fna"), emit: fasta
    tuple val(meta), path("*.clstr")    , emit: clusters
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-c 0.95 -n 10 -d 0 -M 0 -T 0'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cd-hit-est \\
        -i ${fasta} \\
        -o group_${prefix}.dedup.fna \\
        ${args} \\
        -T ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: \$(cd-hit-est -h 2>&1 | grep 'CD-HIT version' | sed 's/.*CD-HIT version //g' | tr -d '\\t')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dedup.fna
    touch ${prefix}.dedup.fna.clstr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: \$(cd-hit-est -h 2>&1 | grep 'CD-HIT version' | sed 's/.*CD-HIT version //g' | tr -d '\\t')
    END_VERSIONS
    """
}
