process EXTRACT_ARG_SEQUENCES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.2--h9ee0642_0':
        'biocontainers/seqkit:2.8.2--h9ee0642_0' }"

    input:
    tuple val(meta), path(amrfinder_tsv), path(annotation_fna)

    output:
    tuple val(meta), path("*.arg_sequences.fna"), emit: arg_sequences
    tuple val(meta), path("*.arg_summary.tsv")  , emit: summary
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # AMRFinderPlus protein mode with Pyrodigal output columns:
    # Column 1: Protein id (gene ID from Pyrodigal, e.g., k141_93725_5)
    # Column 2: Element symbol (gene name, e.g., vanR)
    # Column 7: Class (e.g., GLYCOPEPTIDE)
    # Column 8: Subclass (e.g., VANCOMYCIN)
    
    # Extract gene IDs from AMRFinderPlus output (column 1)
    tail -n +2 ${amrfinder_tsv} | cut -f1 > ${prefix}_arg_ids.txt

    # Check if any ARGs were found
    if [ ! -s ${prefix}_arg_ids.txt ]; then
        echo "No ARGs found in AMRFinderPlus results, creating empty output files"
        touch ${prefix}.arg_sequences.fna
        echo -e "Gene_ID\\tGene_Symbol\\tClass\\tSubclass" > ${prefix}.arg_summary.tsv
    else
        # Create summary TSV with gene details
        echo -e "Gene_ID\\tGene_Symbol\\tClass\\tSubclass" > ${prefix}.arg_summary.tsv
        tail -n +2 ${amrfinder_tsv} | awk -F'\\t' '{print \$1"\\t"\$2"\\t"\$7"\\t"\$8}' >> ${prefix}.arg_summary.tsv

        # Extract NUCLEOTIDE sequences from FNA using gene IDs
        seqkit grep -f ${prefix}_arg_ids.txt ${annotation_fna} > ${prefix}_temp.fna

        # Debug: show what we're matching
        echo "=== First 3 lines of summary TSV ===" >&2
        head -n 3 ${prefix}.arg_summary.tsv >&2
        
        echo "=== First 3 gene IDs from TSV ===" >&2
        head -n 3 ${prefix}_arg_ids.txt >&2
        
        echo "=== First header from temp FNA ===" >&2
        head -n 1 ${prefix}_temp.fna >&2
        
        echo "=== First gene ID extracted from FNA header ===" >&2
        head -n 1 ${prefix}_temp.fna | sed 's/^>//; s/ .*//' >&2

        # Add ARG annotations to FASTA headers
        awk -F'\\t' '
        NR==FNR {
            if (NR > 1) {
                # Read summary: Gene_ID, Gene_Symbol, Class, Subclass
                gene[\$1] = \$2;
                class[\$1] = \$3;
                subclass[\$1] = \$4;
                print "DEBUG: Loaded gene_id=" \$1 " gene_symbol=" \$2 > "/dev/stderr";
            }
            next;
        }
        /^>/ {
            # Extract gene ID (first word after >) - simpler approach for BusyBox awk
            line = \$0;
            sub(/^>/, "", line);
            sub(/ .*/, "", line);
            gene_id = line;
            print "DEBUG: Looking up gene_id=" gene_id > "/dev/stderr";
            
            if (gene_id in gene) {
                header = ">" gene_id " | " gene[gene_id] " | " class[gene_id];
                if (subclass[gene_id] != "" && subclass[gene_id] != "NA") {
                    header = header " | " subclass[gene_id];
                }
                print "DEBUG: MATCHED - printing new header" > "/dev/stderr";
                print header;
            } else {
                print "DEBUG: NOT FOUND - keeping original" > "/dev/stderr";
                print \$0;
            }
            next;
        }
        {print}
        ' ${prefix}.arg_summary.tsv ${prefix}_temp.fna > ${prefix}.arg_sequences.fna

        rm ${prefix}_temp.fna
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.arg_sequences.fna
    touch ${prefix}.arg_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//g')
    END_VERSIONS
    """
}
