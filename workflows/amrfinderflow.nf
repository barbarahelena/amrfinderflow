/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ANNOTATION                } from '../subworkflows/local/annotation'
include { PROTEIN_ANNOTATION        } from '../subworkflows/local/protein_annotation'
include { ARG                       } from '../subworkflows/local/arg'
include { TAXA_CLASS                } from '../subworkflows/local/taxa_class'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_INPUT_PREP     } from '../modules/nf-core/gunzip/main'
include { SEQKIT_SEQ as SEQKIT_SEQ_LENGTH } from '../modules/nf-core/seqkit/seq/main'
include { SEQKIT_SEQ as SEQKIT_SEQ_FILTER } from '../modules/nf-core/seqkit/seq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AMRFINDERFLOW {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_reads        // channel: reads read in from --fastqs, for mapping back to ARG catalog if required

    main:

    ch_versions = Channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INPUT PREPARATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // Some tools require uncompressed input
    ch_input_prep = ch_samplesheet
        .map { meta, fasta, faa, gbk -> [meta + [category: 'all'], [fasta, faa, gbk]] }
        .transpose()
        .branch {
            compressed: it[1].toString().endsWith('.gz')
            uncompressed: it[1]
        }

    GUNZIP_INPUT_PREP(ch_input_prep.compressed)
    ch_versions = ch_versions.mix(GUNZIP_INPUT_PREP.out.versions)

    // Merge all the already uncompressed and newly compressed FASTAs here into
    // a single input channel for downstream
    ch_intermediate_input = GUNZIP_INPUT_PREP.out.gunzip
        .mix(ch_input_prep.uncompressed)
        .groupTuple()
        .map { meta, files ->
            def fasta_found = files.find { it.toString().tokenize('.').last().matches('fasta|fas|fna|fa') }
            def faa_found = files.find { it.toString().endsWith('.faa') }
            def gbk_found = files.find { it.toString().tokenize('.').last().matches('gbk|gbff') }
            def fasta = fasta_found != null ? fasta_found : []
            def faa = faa_found != null ? faa_found : []
            def gbk = gbk_found != null ? gbk_found : []

            [meta, fasta, faa, gbk]
        }
        .branch { meta, fasta, faa, gbk ->
            preannotated: gbk != []
            fastas: true
        }

    ch_input_for_annotation = ch_intermediate_input.fastas.map { meta, fasta, protein, gbk -> [meta, fasta] }

    /*
        ANNOTATION
    */

    // Annotation is required for ARG screening with AMRFinderPlus
    if (params.run_arg_screening || params.run_amp_screening || params.run_protein_annotation) {
        ANNOTATION(ch_input_for_annotation)
        ch_versions = ch_versions.mix(ANNOTATION.out.versions)

        ch_new_annotation = ch_input_for_annotation
            .join(ANNOTATION.out.faa)
            .join(ANNOTATION.out.fna)
            .join(ANNOTATION.out.gbk)
    }
    else {
        ch_new_annotation = ch_intermediate_input.fastas
    }

    // Mix back the preannotated samples with the newly annotated ones
    ch_prepped_input = ch_new_annotation
        .mix(ch_intermediate_input.preannotated)
        .multiMap { meta, fasta, faa, fna, gbk ->
            fastas: [meta, fasta]
            faas: [meta, faa]
            fnas: [meta, fna]
            gbks: [meta, gbk]
        }

    /*
        TAXONOMIC CLASSIFICATION
    */

    // The final subworkflow reports need taxonomic classification.
    // This can be either on NT or AA level depending on annotation.
    // TODO: Only NT at the moment. AA tax. classification will be added only when its PR is merged.
    if (params.run_taxa_classification) {
        TAXA_CLASS(ch_prepped_input.fastas)
        ch_versions = ch_versions.mix(TAXA_CLASS.out.versions)
        ch_taxonomy_tsv = TAXA_CLASS.out.sample_taxonomy
    }
    else {

        ch_mmseqs_db = Channel.empty()
        ch_taxonomy_querydb = Channel.empty()
        ch_taxonomy_querydb_taxdb = Channel.empty()
        ch_taxonomy_tsv = Channel.empty()
    }

    /*
        PROTEIN ANNOTATION
    */
    if (params.run_protein_annotation) {
        def filtered_faas = ch_prepped_input.faas.filter { meta, file ->
            if (file != [] && file.isEmpty()) {
                log.warn("[barbarahelena/amrfinderflow] Annotation of the following sample produced an empty FAA file. InterProScan classification of the CDS requiring this file will not be executed: ${meta.id}")
            }
            !file.isEmpty()
        }

        SEQKIT_SEQ_FILTER(filtered_faas)
        ch_versions = ch_versions.mix(SEQKIT_SEQ_FILTER.out.versions)
        ch_input_for_protein_annotation =  SEQKIT_SEQ_FILTER.out.fastx

        PROTEIN_ANNOTATION ( ch_input_for_protein_annotation )
        ch_versions = ch_versions.mix(PROTEIN_ANNOTATION.out.versions)

        ch_interproscan_tsv = PROTEIN_ANNOTATION.out.tsv.map { meta, file ->
            if (file == [] || file.isEmpty()) {
                log.warn("[barbarahelena/amrfinderflow] Protein annotation with InterProScan produced an empty TSV file. No protein annotation will be added for sample ${meta.id}.")
                [meta, []]
            } else {
                [meta, file]
            }
        }
    } else {
        ch_interproscan_tsv = ch_prepped_input.faas.map { meta, _ ->
            [meta, []]
        }
    }


    /*
        SCREENING
    */

    /*
        ARGs - Simplified to only AMRFinderPlus
    */
    if (params.run_arg_screening && !params.run_taxa_classification) {
        ARG(
            ch_prepped_input.fastas,
            ch_prepped_input.faas,
            ch_prepped_input.fnas,
            ch_taxonomy_tsv,
            ch_reads
        )
        ch_versions = ch_versions.mix(ARG.out.versions)
    }
    else if (params.run_arg_screening && params.run_taxa_classification) {
        ARG(
            ch_prepped_input.fastas,
            ch_prepped_input.faas,
            ch_prepped_input.fnas,
            ch_taxonomy_tsv.filter { meta, file ->
                if (file.isEmpty()) {
                    log.warn("[barbarahelena/amrfinderflow] Taxonomy classification of the following sample produced an empty TSV file. Taxonomy merging will not be executed: ${meta.id}")
                }
                !file.isEmpty()
            },
            ch_reads
        )
        ch_versions = ch_versions.mix(ARG.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'amrfinderflow_software_versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
