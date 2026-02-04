//
// Subworkflow with functionality specific to the barbarahelena/amrfinderflow pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification          } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()
    
    //
    // Create channel from input file provided through params.input
    //
    Channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .set { ch_samplesheet }
    
    //
    // Create channel from fastqs file if provided (for read mapping back to ARG catalog)
    //
    if (params.fastqs) {
        Channel.fromList(samplesheetToList(params.fastqs, "${projectDir}/assets/schema_fastqs.json"))
            .map { row -> validateFastqsSamplesheet(row) }
            .set { ch_reads }
    } else {
        ch_reads = Channel.empty()
    }


    emit:
    samplesheet = ch_samplesheet
    reads = ch_reads
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
}

//
// Validate channels from fastqs sheet
//
def validateFastqsSamplesheet(input) {
    // Input format from samplesheetToList: [meta, fastq_1, fastq_2]
    // where meta contains: [id: sample, group: group, single_end: false]
    def meta = input[0]
    def fastq_1 = input[1]
    def fastq_2 = input[2]

    // Check that fastq files exist
    if (fastq_1 && !fastq_1.exists()) {
        error("Please check fastqs samplesheet -> FastQ file does not exist: ${fastq_1}")
    }
    
    if (fastq_2 && !fastq_2.exists()) {
        error("Please check fastqs samplesheet -> FastQ file does not exist: ${fastq_2}")
    }

    return [meta, [fastq_1, fastq_2]]
}

//
// Generate methods description
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    def preprocessing_text = "The pipeline used the following tools: preprocessing included SeqKit2 (Shen et al. 2024)."

    def annotation_text = [
        "Annotation was carried out with:",
        params.annotation_tool == 'prodigal' ? "Prodigal (Hyatt et al. 2010)." : "",
        params.annotation_tool == 'pyrodigal' ? "Pyrodigal (Larralde 2022)." : "",
        params.annotation_tool == 'bakta' ? "BAKTA (Schwengers et al. 2021)." : "",
        params.annotation_tool == 'prokka' ? "PROKKA (Seemann 2014)." : "",
    ].join(' ').trim()

    def amp_text = [
        "The following antimicrobial peptide screening tools were used:",
        !params.amp_skip_amplify ? "AMPlify (Li et al. 2022)," : "",
        !params.amp_skip_macrel ? "Macrel (Santos-Júnior et al. 2020)," : "",
        !params.amp_skip_ampir ? "ampir (Fingerhut et al. 2021)," : "",
        params.amp_run_hmmsearch ? "HMMER (Eddy 2011)," : "",
        ". The output from the antimicrobial peptide screening tools were standardised and summarised with AMPcombi (Ibrahim and Perelo 2023).",
    ].join(' ').trim().replaceAll(', .', ".")

    def arg_text = [
        "The following antimicrobial resistance gene screening tools were used:",
        !params.arg_skip_fargene ? "fARGene (Berglund et al. 2019)," : "",
        !params.arg_skip_rgi ? "RGI (Alcock et al. 2020)," : "",
        !params.arg_skip_amrfinderplus ? "AMRfinderplus (Feldgarden et al. 2021)," : "",
        !params.arg_skip_deeparg ? "deepARG (Arango-Argoty 2018)," : "",
        !params.arg_skip_abricate ? "ABRicate (Seemann 2020)," : "",
        !params.arg_skip_argnorm ? ". The outputs from ARG screening tools were normalized to the antibiotic resistance ontology using argNorm (Ugarcina Perovic et al. 2025)," : "",
        ". The output from the antimicrobial resistance gene screening tools were standardised and summarised with hAMRonization (Maguire et al. 2023).",
    ].join(' ').trim().replaceAll(', +.', ".")

    def bgc_text = [
        "The following biosynthetic gene cluster screening tools were used:",
        !params.bgc_skip_antismash ? "antiSMASH (Blin et al. 2021)," : "",
        !params.bgc_skip_deepbgc ? "deepBGC (Hannigan et al. 2019)," : "",
        !params.bgc_skip_gecco ? "GECCO (Carroll et al. 2021)," : "",
        params.bgc_run_hmmsearch ? "HMMER (Eddy 2011)," : "",
        ". The output from the biosynthetic gene cluster screening tools were standardised and summarised with comBGC (Frangenberg et al. 2023).",
    ].join(' ').replaceAll(', +.', ".").trim()

    def postprocessing_text = "Run statistics were reported using MultiQC (Ewels et al. 2016)."

    def citation_text = [
        preprocessing_text,
        annotation_text,
        params.run_amp_screening ? amp_text : "",
        params.run_arg_screening ? arg_text : "",
        params.run_bgc_screening ? bgc_text : "",
        postprocessing_text,
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? '<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    def preprocessing_text = '<li>Shen, W., Sipos, B., & Zhao, L. (2024). SeqKit2: A Swiss army knife for sequence and alignment processing. iMeta, e191. <a href="https://doi.org/10.1002/imt2.191">https://doi.org/10.1002/imt2.191</a></li>'

    def annotation_text = [
        params.annotation_tool == 'prodigal' ? '<li>Hyatt, D., Chen, G. L., Locascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics, 11, 119. DOI: <a href="https://doi.org/10.1186/1471-2105-11-119">10.1186/1471-2105-11-119</a></li>' : "",
        params.annotation_tool == 'pyrodigal' ? '<li>Larralde, M. (2022). Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. Journal of Open Source Software, 7(72), 4296. DOI: <a href="https://doi.org/10.21105/joss.04296">10.21105/joss.04296</a></li>' : "",
        params.annotation_tool == 'bakta' ? '<li>Schwengers, O., Jelonek, L., Dieckmann, M. A., Beyvers, S., Blom, J., & Goesmann, A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11). DOI: <a href="https://doi.org/10.1099/mgen.0.000685">10.1099/mgen.0.000685</a></li>' : "",
        params.annotation_tool == 'prokka' ? '<li>Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics (Oxford, England), 30(14), 2068–2069. DOI: <a href="https://doi.org/10.1093/bioinformatics/btu153">10.1093/bioinformatics/btu153</a></li>' : "",
    ].join(' ').trim()
    def arg_text = [
        !params.arg_skip_fargene ? '<li>Berglund, F., Österlund, T., Boulund, F., Marathe, N. P., Larsson, D., & Kristiansson, E. (2019). Identification and reconstruction of novel antibiotic resistance genes from metagenomes. Microbiome, 7(1), 52. DOI: <a href="https://doi.org/10.1186/s40168-019-0670-1">10.1186/s40168-019-0670-1</a></li>' : "",
        !params.arg_skip_rgi ? '<li>Alcock, B. P., Raphenya, A. R., Lau, T., Tsang, K. K., Bouchard, M., Edalatmand, A., Huynh, W., Nguyen, A. V., Cheng, A. A., Liu, S., Min, S. Y., Miroshnichenko, A., Tran, H. K., Werfalli, R. E., Nasir, J. A., Oloni, M., Speicher, D. J., Florescu, A., Singh, B., Faltyn, M., … McArthur, A. G. (2020). CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. Nucleic acids research, 48(D1), D517–D525. DOI: <a href="https://doi.org/10.1093/nar/gkz935">10.1093/nar/gkz935</a></li>' : "",
        !params.arg_skip_amrfinderplus ? '<li>Feldgarden, M., Brover, V., Gonzalez-Escalona, N., Frye, J. G., Haendiges, J., Haft, D. H., Hoffmann, M., Pettengill, J. B., Prasad, A. B., Tillman, G. E., Tyson, G. H., & Klimke, W. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Scientific reports, 11(1), 12728. DOI: <a href="https://doi.org/10.1038/s41598-021-91456-0">10.1038/s41598-021-91456-0</a></li>' : "",
        !params.arg_skip_deeparg ? '<li>Arango-Argoty, G., Garner, E., Pruden, A., Heath, L. S., Vikesland, P., & Zhang, L. (2018). DeepARG: a deep learning approach for predicting antibiotic resistance genes from metagenomic data. Microbiome, 6(1), 23. DOI: <a href="https://doi.org/10.1186/s40168-018-0401-z">10.1186/s40168-018-0401-z</a></li>' : "",
        !params.arg_skip_abricate ? '<li>Seemann, T. (2020). ABRicate. Github <a href="https://github.com/tseemann/abricate">https://github.com/tseemann/abricate</a>.</li>' : "",
        !params.arg_skip_argnorm ? '<li>Ugarcina Perovic, S., Ramji, V., Chong, H., Duan, Y., Maguire, F., Coelho, L. P. (2025). argNorm: normalization of antibiotic resistance gene annotations to the Antibiotic Resistance Ontology (ARO), Bioinformatics, btaf173. DOI: <a href="https://doi.org/10.1093/bioinformatics/btaf173">10.1093/bioinformatics/btaf173</a></li>' : "",
        '<li>Public Health Alliance for Genomic Epidemiology (pha4ge). (2022). Parse multiple Antimicrobial Resistance Analysis Reports into a common data structure. Github. Retrieved October 5, 2022, from <a href="https://github.com/pha4ge/hAMRonization">https://github.com/pha4ge/hAMRonization</a></li>',
    ].join(' ').trim().replaceAll(', +.', ".")

    def postprocessing_text = '<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. <a href="https://doi.org/10.1093/bioinformatics/btw354">https://doi.org/10.1093/bioinformatics/btw354</a></li>'

    // Special as reused in multiple subworkflows, and we don't want to cause duplicates
    def reference_text = [
        preprocessing_text,
        annotation_text,
        arg_text,
        postprocessing_text,
    ].join(' ').trim()

    return reference_text
}