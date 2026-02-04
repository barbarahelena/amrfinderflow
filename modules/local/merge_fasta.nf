process MERGE_FASTA {
    tag "Group $meta.group"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.merged.faa"), emit: fasta
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.group}"
    """
    # Merge all FASTA files
    cat ${fastas} > group_${prefix}.merged.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version | head -n1 | sed 's/cat (GNU coreutils) //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version | head -n1 | sed 's/cat (GNU coreutils) //g')
    END_VERSIONS
    """
}
