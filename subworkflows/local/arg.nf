/*
    Run ARG screening tools - AMRFinderPlus in protein mode
    Uses protein sequences for better sensitivity, extracts nucleotide sequences for read mapping
*/

include { AMRFINDERPLUS_UPDATE                          } from '../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN                             } from '../../modules/nf-core/amrfinderplus/run/main'
include { EXTRACT_ARG_SEQUENCES                         } from '../../modules/local/extract_arg_sequences'
include { MERGE_FASTA                                   } from '../../modules/local/merge_fasta'
include { CDHIT_EST                                     } from '../../modules/local/cdhit_est'
include { BWA_MEM2_INDEX                                } from '../../modules/local/bwa_mem2_index'
include { BWA_MEM2_MAP_ARG                              } from '../../modules/local/bwa_mem2_map_arg'
include { MERGE_ARG_COUNTS as MERGE_ARG_COUNTS_GROUP    } from '../../modules/local/merge_arg_counts'
include { MERGE_ARG_COUNTS as MERGE_ARG_COUNTS_ALL      } from '../../modules/local/merge_arg_counts'

workflow ARG {
    take:
    _fastas         // tuple val(meta), path(contigs) - not used in protein mode
    faa_annotations // tuple val(meta), path(faa) - Pyrodigal protein files for AMRFinderPlus
    fna_annotations // tuple val(meta), path(fna) - Pyrodigal nucleotide sequences for extraction
    _tsvs           // tuple val(meta), path(MMSEQS_CREATETSV.out.tsv) - not used
    reads           // tuple val(meta), path(reads) - for mapping back

    main:
    ch_versions = Channel.empty()

    // AMRfinderplus run in protein mode for better sensitivity
    // Prepare channel for database
    if (params.arg_amrfinderplus_db) {
        ch_amrfinderplus_db = Channel
            .fromPath(params.arg_amrfinderplus_db, checkIfExists: true)
            .first()
    }
    else {
        AMRFINDERPLUS_UPDATE()
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
        ch_amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
    }

    // Add is_proteins flag to metadata for protein mode
    ch_amrfinder_input = faa_annotations.map { meta, faa ->
        def new_meta = meta + [is_proteins: true]
        [new_meta, faa]
    }

    AMRFINDERPLUS_RUN(ch_amrfinder_input, ch_amrfinderplus_db)
    ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

    // Extract ARG nucleotide sequences from FNA annotations
    // Need to remove is_proteins flag from metadata for join to work
    ch_extract_input = AMRFINDERPLUS_RUN.out.report
        .map { meta, tsv ->
            // Remove is_proteins flag to match original annotations metadata
            def clean_meta = meta.findAll { k, v -> k != 'is_proteins' }
            [clean_meta, tsv]
        }
        .join(fna_annotations)
    
    EXTRACT_ARG_SEQUENCES(ch_extract_input)
    ch_versions = ch_versions.mix(EXTRACT_ARG_SEQUENCES.out.versions)

    // Group by group if provided, otherwise by sample
    // meta.group will be used to group samples, or fall back to sample ID
    ch_arg_sequences = EXTRACT_ARG_SEQUENCES.out.arg_sequences
        .map { meta, fasta ->
            // Use group if it exists (even if it's 0), otherwise use sample id
            def group_id = meta.containsKey('group') && meta.group != null ? meta.group : meta.id
            [group_id, meta, fasta]
        }
        .groupTuple(by: 0)
        .map { group_id, metas, fastas ->
            // Merge metadata from all samples in this group
            def all_samples = metas.collect { m -> m.id }
            def merged_meta = [id: group_id, group: group_id, samples: all_samples]
            [merged_meta, fastas]
        }

    // Merge FASTAs per group
    MERGE_FASTA(ch_arg_sequences)
    ch_versions = ch_versions.mix(MERGE_FASTA.out.versions)

    // Deduplicate with CD-HIT-EST
    CDHIT_EST(MERGE_FASTA.out.fasta)
    ch_versions = ch_versions.mix(CDHIT_EST.out.versions)

    // Index deduplicated ARG catalogs
    BWA_MEM2_INDEX(CDHIT_EST.out.fasta)
    ch_versions = ch_versions.mix(BWA_MEM2_INDEX.out.versions)

    // Map reads back to ARG catalog
    // Match reads to catalogs by group field
    // Group FNA files per group (same logic as ARG sequences)
    
    ch_fna_grouped = fna_annotations
        .map { meta, fna ->
            def group_id = meta.containsKey('group') && meta.group != null ? meta.group : meta.id
            [group_id, fna]
        }
        .groupTuple(by: 0)
    
    // First join catalog index with grouped FNA files by group
    ch_catalog_with_fna = BWA_MEM2_INDEX.out.index
        .map { meta, fasta, index_files ->
            [meta.group, fasta, index_files]
        }
        .join(ch_fna_grouped, by: 0)
        .map { group, fasta, index_files, fna_files ->
            // Return: [group, fasta, index_files, fna_files]
            [group, fasta, index_files, fna_files]
        }
    
    // Now combine with each individual read sample in that group
    ch_catalog_with_reads = ch_catalog_with_fna
        .combine(reads.map { meta, files -> [meta.group, meta, files] }, by: 0)
        .map { _group, fasta, index_files, fna_files, read_meta, read_files ->
            // Return format: [read_meta, reads, catalog_fasta, index_files, fna_files]
            // Each sample in the group gets the same catalog and FNA files
            [read_meta, read_files, fasta, index_files, fna_files]
        }

    BWA_MEM2_MAP_ARG(ch_catalog_with_reads)
    ch_versions = ch_versions.mix(BWA_MEM2_MAP_ARG.out.versions)

    // Merge arg_counts.tsv files per group
    ch_counts_grouped = BWA_MEM2_MAP_ARG.out.counts
        .map { meta, counts ->
            def group_id = meta.containsKey('group') && meta.group != null ? meta.group : meta.id
            [group_id, meta, counts]
        }
        .groupTuple(by: 0)
        .map { group_id, _metas, counts_files ->
            def merged_meta = [id: group_id, group: group_id]
            [merged_meta, counts_files]
        }

    MERGE_ARG_COUNTS_GROUP(ch_counts_grouped)
    ch_versions = ch_versions.mix(MERGE_ARG_COUNTS_GROUP.out.versions)

    // Also create a population-wide merge of all samples
    ch_counts_all = BWA_MEM2_MAP_ARG.out.counts
        .map { _meta, counts -> counts }
        .collect()
        .map { all_counts ->
            def population_meta = [id: 'all_samples', group: 'population']
            [population_meta, all_counts]
        }

    MERGE_ARG_COUNTS_ALL(ch_counts_all)
    ch_versions = ch_versions.mix(MERGE_ARG_COUNTS_ALL.out.versions)

    // Reporting - Skip hamronization in protein mode (incompatible output format)
    // Hamronization expects nucleotide mode output with "Contig id" column
    // In protein mode, AMRFinderPlus output format differs and hamronization fails
    // Direct AMRFinderPlus output is still available for downstream use

    emit:
    versions                = ch_versions
    arg_bam                 = BWA_MEM2_MAP_ARG.out.bam
    arg_idxstats            = BWA_MEM2_MAP_ARG.out.idxstats
    arg_counts              = BWA_MEM2_MAP_ARG.out.counts
    arg_counts_merged_group = MERGE_ARG_COUNTS_GROUP.out.merged_counts
    arg_counts_merged_all   = MERGE_ARG_COUNTS_ALL.out.merged_counts
    arg_catalog             = CDHIT_EST.out.fasta
}
