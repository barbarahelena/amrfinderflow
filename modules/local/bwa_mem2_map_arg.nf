process BWA_MEM2_MAP_ARG {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bwa-mem2=2.3 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bwa-mem2_samtools:4cdd9550c1e63945':
        'community.wave.seqera.io/library/bwa-mem2_samtools:47cbe92091673673' }"

    input:
    tuple val(meta), path(reads), path(catalog), path(index), path(fna)

    output:
    tuple val(meta), path("*.bam")     , emit: bam
    tuple val(meta), path("*.bam.bai") , emit: bai
    tuple val(meta), path("*.idxstats.tsv"), emit: idxstats
    tuple val(meta), path("*.arg_counts.tsv"), emit: counts
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_args = task.ext.args2 ?: '-b -F 4 -q 10'
    """
    # Map reads to ARG catalog
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        ${args} \\
        ${catalog} \\
        ${reads[0]} \\
        ${reads[1]} \\
        | samtools view ${samtools_args} \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.bam

    # Index BAM
    samtools index ${prefix}.bam

    # Get total number of reads in the sample (both mapped and unmapped to catalog)
    TOTAL_READS=\$(samtools view -c -F 0x900 ${reads[0]})
    
    # Generate mapping statistics
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats.tsv

    # Count total genes predicted by Pyrodigal across all FNA files in the group
    TOTAL_GENES=0
    for fna_file in ${fna}; do
        GENE_COUNT=\$(grep -c "^>" \$fna_file || echo "0")
        TOTAL_GENES=\$((TOTAL_GENES + GENE_COUNT))
    done

    # Create tidy TSV with gene metadata from FASTA headers
    echo -e "Sample\\tGene_ID\\tGene_Symbol\\tClass\\tSubclass\\tSequence_Length\\tMapped_Reads\\tUnmapped_Reads\\tTotal_Genes\\tTotal_Reads\\tPrevalence\\tRPK\\tRPKM\\tCoverage" > ${prefix}.arg_counts.tsv
    
    # Parse idxstats and extract metadata from catalog headers
    awk -v total_genes=\$TOTAL_GENES -v total_reads=\$TOTAL_READS 'BEGIN {FS="\\t"; OFS="\\t"}
    NR==FNR {
        if (\$0 ~ /^>/) {
            # Remove > and get full header
            header = \$0
            sub(/^>/, "", header)
            
            # Split by pipe to get parts: gene_id | gene_symbol | class | subclass
            n = split(header, parts, " \\\\| ")
            
            # Extract gene_id (first part, remove trailing spaces)
            gene_id = parts[1]
            gsub(/^ +| +\$/, "", gene_id)
            
            # Extract remaining fields with defaults
            gene_symbol = (n >= 2) ? parts[2] : "NA"
            gsub(/^ +| +\$/, "", gene_symbol)
            
            classval = (n >= 3) ? parts[3] : "NA"
            gsub(/^ +| +\$/, "", classval)
            
            subclass = (n >= 4) ? parts[4] : "NA"
            gsub(/^ +| +\$/, "", subclass)
            
            # Store metadata
            meta[gene_id] = gene_symbol "\\t" classval "\\t" subclass
        }
        next
    }
    \$1 != "*" {
        # Process idxstats
        gene_id = \$1
        len = \$2
        mapped = \$3
        unmapped = \$4
        
        # Prevalence: present (1) if any reads mapped, absent (0) otherwise
        prevalence = (mapped > 0) ? 1 : 0
        
        # Calculate RPK: (reads * 1000) / gene_length
        # This is not normalized for sequencing depth - useful for within-sample comparisons
        if (len > 0) {
            rpk = (mapped * 1000) / len
        } else {
            rpk = 0
        }
        
        # Calculate RPKM based on total sample reads: (reads * 10^9) / (gene_length * total_sample_reads)
        # This normalizes for both sequencing depth and gene length
        if (len > 0 && total_reads > 0) {
            rpkm = (mapped * 1000000000) / (len * total_reads)
        } else {
            rpkm = 0
        }
        
        # Calculate coverage: (reads * read_length) / gene_length
        # Assuming 150bp reads (standard Illumina)
        if (len > 0) {
            coverage = (mapped * 150) / len
        } else {
            coverage = 0
        }
        
        if (gene_id in meta) {
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%d\\t%.2f\\t%.2f\\t%.2f\\n", 
                "${meta.id}", gene_id, meta[gene_id], len, mapped, unmapped, 
                total_genes, total_reads, prevalence, rpk, rpkm, coverage
        } else {
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%d\\t%.2f\\t%.2f\\t%.2f\\n", 
                "${meta.id}", gene_id, "NA", "NA", "NA", len, mapped, unmapped, 
                total_genes, total_reads, prevalence, rpk, rpkm, coverage
        }
    }' ${catalog} ${prefix}.idxstats.tsv >> ${prefix}.arg_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version | head -n1 | sed 's/bwa-mem2 //g')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}.idxstats.tsv
    touch ${prefix}.arg_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version | head -n1 | sed 's/bwa-mem2 //g')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //g')
    END_VERSIONS
    """
}
