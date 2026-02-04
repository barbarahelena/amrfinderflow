process MERGE_ARG_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(counts, stageAs: 'counts_*.tsv')

    output:
    tuple val(meta), path("*.merged_arg_counts.tsv"), emit: merged_counts
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Get all count files (handles both counts_.tsv and counts_N.tsv patterns)
    COUNT_FILES=\$(ls counts*.tsv 2>/dev/null)
    
    if [ -z "\$COUNT_FILES" ]; then
        echo "ERROR: No count files found"
        ls -la
        exit 1
    fi
    
    # Get header from first file
    FIRST_FILE=\$(echo "\$COUNT_FILES" | head -n 1)
    head -n 1 "\$FIRST_FILE" > ${prefix}.merged_arg_counts.tsv
    
    # Append all data (skip headers) from all files
    for file in \$COUNT_FILES; do
        tail -n +2 "\$file" >> ${prefix}.merged_arg_counts.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(head --version | head -n1 | sed 's/head (GNU coreutils) //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged_arg_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(head --version | head -n1 | sed 's/head (GNU coreutils) //')
    END_VERSIONS
    """
}
