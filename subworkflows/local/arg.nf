/*
    Run ARG screening tools - Simplified to AMRFinderPlus + argNorm only
*/

include { AMRFINDERPLUS_UPDATE             } from '../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN                } from '../../modules/nf-core/amrfinderplus/run/main'
include { TABIX_BGZIP as ARG_TABIX_BGZIP   } from '../../modules/nf-core/tabix/bgzip/main'
include { MERGE_TAXONOMY_HAMRONIZATION     } from '../../modules/local/merge_taxonomy_hamronization'
include { HAMRONIZATION_SUMMARIZE          } from '../../modules/nf-core/hamronization/summarize/main'
include { HAMRONIZATION_AMRFINDERPLUS      } from '../../modules/nf-core/hamronization/amrfinderplus/main'
include { ARGNORM as ARGNORM_AMRFINDERPLUS } from '../../modules/nf-core/argnorm/main'

workflow ARG {
    take:
    fastas      // tuple val(meta), path(contigs)
    annotations
    tsvs        // tuple val(meta), path(MMSEQS_CREATETSV.out.tsv)

    main:
    ch_versions = Channel.empty()

    // Prepare HAMRONIZATION reporting channel
    ch_input_to_hamronization_summarize = Channel.empty()

    // AMRfinderplus run
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

    AMRFINDERPLUS_RUN(fastas, ch_amrfinderplus_db)
    ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

    // Reporting
    HAMRONIZATION_AMRFINDERPLUS(AMRFINDERPLUS_RUN.out.report, 'tsv', AMRFINDERPLUS_RUN.out.tool_version, AMRFINDERPLUS_RUN.out.db_version)
    ch_versions = ch_versions.mix(HAMRONIZATION_AMRFINDERPLUS.out.versions)
    ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_AMRFINDERPLUS.out.tsv)

    if (!params.arg_skip_argnorm) {
        ch_input_to_argnorm_amrfinderplus = HAMRONIZATION_AMRFINDERPLUS.out.tsv.filter { meta, file -> !file.isEmpty() }
        ARGNORM_AMRFINDERPLUS(ch_input_to_argnorm_amrfinderplus, 'amrfinderplus', 'ncbi')
        ch_versions = ch_versions.mix(ARGNORM_AMRFINDERPLUS.out.versions)
    }

    ch_input_to_hamronization_summarize
        .map {
            it[1]
        }
        .collect()
        .set { ch_input_for_hamronization_summarize }

    HAMRONIZATION_SUMMARIZE(ch_input_for_hamronization_summarize, params.arg_hamronization_summarizeformat)
    ch_versions = ch_versions.mix(HAMRONIZATION_SUMMARIZE.out.versions)

    // MERGE_TAXONOMY
    if (params.run_taxa_classification) {

        ch_mmseqs_taxonomy_list = tsvs.map { it[1] }.collect()
        MERGE_TAXONOMY_HAMRONIZATION(HAMRONIZATION_SUMMARIZE.out.tsv, ch_mmseqs_taxonomy_list)
        ch_versions = ch_versions.mix(MERGE_TAXONOMY_HAMRONIZATION.out.versions)

        ch_tabix_input = Channel
            .of(['id': 'hamronization_combined_report'])
            .combine(MERGE_TAXONOMY_HAMRONIZATION.out.tsv)

        ARG_TABIX_BGZIP(ch_tabix_input)
        ch_versions = ch_versions.mix(ARG_TABIX_BGZIP.out.versions)
    }

    emit:
    versions = ch_versions
}
