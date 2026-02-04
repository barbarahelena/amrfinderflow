# AMRfinderflow: Output

## Introduction

AMRfinderflow is a simplified pipeline focused on **antimicrobial resistance gene (ARG) detection** using:

- [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder) – antimicrobial resistance gene detection, using NCBI's curated Reference Gene Database and curated collection of Hidden Markov Models
- [hAMRonization](https://github.com/pha4ge/hAMRonization) – standardized reporting of AMR results
- [argNorm](https://github.com/BigDataBiology/argNorm) – normalization of ARG annotations to the [ARO](https://obofoundry.org/ontology/aro.html) (Antibiotic Resistance Ontology)

The pipeline workflow includes:
1. Sequence quality control with [SeqKit](https://bioinf.shenwei.me/seqkit/)
2. Optional taxonomic classification with [MMseqs2](https://github.com/soedinglab/MMseqs2)
3. Gene prediction using [Pyrodigal](https://github.com/althonos/pyrodigal), [Prodigal](https://github.com/hyattpd/Prodigal), [Prokka](https://github.com/tseemann/prokka), or [Bakta](https://github.com/oschwengers/bakta)
4. Optional protein domain annotation with [InterProScan](https://github.com/ebi-pf-team/interproscan)
5. ARG detection with AMRFinderPlus, standardization with hAMRonization, and normalization with argNorm

We recommend first looking at the summary report ([hAMRonization](#hamronization)) to get an overview of ARG hits found, then exploring the specific output directories for detailed information about each result. The tool-specific output directories also include the output from the annotation steps if the `--save_annotations` flag was set. Additionally, taxonomic classifications from [MMseqs2](https://github.com/soedinglab/MMseqs2) are saved if the `--taxa_classification_mmseqs_db_savetmp` and `--taxa_classification_mmseqs_taxonomy_savetmp` flags are set.

All downloaded databases (from [MMseqs2](https://github.com/soedinglab/MMseqs2), [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder), and/or [Bakta](https://github.com/oschwengers/bakta)) are saved into the output directory `<outdir>/databases/` if the `--save_db` flag was set.

The directories listed below will be created in the results directory (specified by the `--outdir` flag) after the pipeline has finished. All paths are relative to this top-level output directory. The default directory structure of AMRfinderflow is:

```tree
results/
├── annotation/
|   ├── bakta/
|   ├── prodigal/
|   ├── prokka/
|   └── pyrodigal/
├── arg/
|   ├── amrfinderplus/
|   ├── argnorm/
|   ├── bwa_mapping/
|   ├── deduplicated_catalog/
|   ├── extracted_sequences/
|   ├── hamronization/
|   ├── merged_counts/
|   └── merged_sequences/
├── databases/
├── pipeline_info/
├── protein_annotation/
|   └── interproscan/
├── qc/
|   └── seqkit/
├── reports/
|   └── hamronization_summarize/
└── taxonomic_classification/
    └── mmseqs_createtsv/
work/
```

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes prokaryotic sequence data through the following steps:

Input contig QC with:

- [SeqKit](https://bioinf.shenwei.me/seqkit/) (default) - for separating into long- and short- categories

Taxonomy classification of nucleotide sequences with:

- [MMseqs2](https://github.com/soedinglab/MMseqs2) (default) - for contig taxonomic classification using 2bLCA.

ORF prediction and annotation with any of:

- [Pyrodigal](#pyrodigal) (default) – for open reading frame prediction.
- [Prodigal](#prodigal) – for open reading frame prediction.
- [Prokka](#prokka) – open reading frame prediction and functional protein annotation.
- [Bakta](#bakta) – open reading frame prediction and functional protein annotation.

CDS domain annotation:

- [InterProScan](#interproscan) (optional) – for open reading frame protein and domain predictions.

Antimicrobial Resistance Genes (ARGs):

- [AMRFinderPlus](#amrfinderplus) – antimicrobial resistance gene detection, using NCBI's curated Reference Gene Database and curated collection of Hidden Markov Models.
- [hAMRonization](#hamronization) – standardized summary of antimicrobial resistance gene output.
- [argNorm](#argnorm) - Normalize ARG annotations from AMRFinderPlus to the ARO (Antibiotic Resistance Ontology).

Output Summaries:

- [hAMRonization](#hamronization) – summary of antimicrobial resistance gene output.
- [argNorm](#argnorm) - Normalized ARG annotations to the ARO.
- [Pipeline information](#pipeline-information) – report metrics generated during the workflow execution.

## Tool details

### Taxonomic classification tools

[MMseqs2](#mmseqs2)

#### MMseqs2

<details markdown="1">
<summary>Output files</summary>

- `taxonomic_classification/mmseqs2_createtsv/`
  - `<samplename>/`:
    - `*.tsv`: tab-separated table containing the taxonomic lineage of every contig. When a contig cannot be classified according to the database, it is assigned in the 'lineage' column as 'no rank | unclassified'.
- `reports/<workflow>/<workflow>_complete_summary_taxonomy.tsv.gz`: tab-separated table containing the concatenated results from the <workflow> summary tables along with the taxonomic classification if the parameter `--run_taxa_classification` is called.
</details>

[MMseqs2](https://github.com/soedinglab/MMseqs2) classifies the taxonomic lineage of contigs based on the last common ancestor. The inferred taxonomic lineages are included in the final workflow summaries to annotate the potential source bacteria of the identified ARGs.

### Annotation tools

[Pyrodigal](#pyrodigal), [Prodigal](#prodigal), [Prokka](#prokka), [Bakta](#bakta)

#### Prodigal

<details markdown="1">
<summary>Output files</summary>

- `prodigal/`
  - `<samplename>/`:
    - `*.fna`: nucleotide FASTA file of the input contig sequences
    - `*.faa`: protein FASTA file of the translated CDS sequences
    - `*.gbk`: annotation in GBK format, containing both sequences and annotations

> Descriptions taken from the [Prodigal documentation](https://github.com/hyattpd/prodigal/wiki)

</details>

[Prodigal](https://github.com/hyattpd/Prodigal) annotates whole (meta-)genomes by identifying ORFs in a set of genomic DNA sequences. The output is used by some of the functional screening tools.

#### Pyrodigal

<details markdown="1">
<summary>Output files</summary>

- `pyrodigal/`
  - `<samplename>/`:
    - `*.gbk`: annotation in GBK format, containing both sequences and annotations
    - `*.fna`: nucleotide FASTA file of the annotated CDS sequences
    - `*.faa`: protein FASTA file of the translated CDS sequences

> Descriptions taken from the [Pyrodigal documentation](https://pyrodigal.readthedocs.io/)

</details>

[Pyrodigal](https://github.com/althonos/pyrodigal) annotates whole (meta-)genomes by identifying ORFs in a set of genomic DNA sequences. It produces the same results as [Prodigal](#prodigal) while being more resource-optimized, thus faster. Unlike Prodigal, Pyrodigal cannot produce output in GenBank format. The output is used by some of the functional screening tools.

#### Prokka

<details markdown="1">
<summary>Output files</summary>

- `prokka/`
  - `<samplename>/`
    - `*.gff`: annotation in GFF3 format, containing both sequences and annotations
    - `*.gbk`: standard Genbank file derived from the master .gff
    - `*.fna`: nucleotide FASTA file of the input contig sequences
    - `*.faa`: protein FASTA file of the translated CDS sequences
    - `*.ffn`: nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
    - `*.sqn`: an ASN1 format "Sequin" file for submission to Genbank
    - `*.fsa`: nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file
    - `*.tbl`: feature Table file, used by "tbl2asn" to create the .sqn file
    - `*.err`: unacceptable annotations - the NCBI discrepancy report
    - `*.log`: logging output that Prokka produced during its run
    - `*.txt`: statistics relating to the annotated features found
      - `*.tsv`: tab-separated file of all features

> Descriptions directly from the [Prokka documentation](https://github.com/tseemann/prokka#output-files)

</details>

[Prokka](https://github.com/tseemann/prokka) performs whole genome annotation to identify features of interest in a set of (meta-)genomic DNA sequences. The output is used by some of the functional screening tools.

#### Bakta

<details markdown="1">
<summary>Output files</summary>

- `bakta/`
  - `<samplename>`
    - `<samplename>.gff3`: annotations & sequences in GFF3 format
    - `<samplename>.gbff`: annotations & sequences in (multi) GenBank format
    - `<samplename>.ffn`: feature nucleotide sequences as FASTA
    - `<samplename>.fna`: replicon/contig DNA sequences as FASTA
    - `<samplename>.embl`: annotations & sequences in (multi) EMBL format
    - `<samplename>.faa`: CDS/sORF amino acid sequences as FASTA
    - `<samplename>_hypothetical.faa`: further information on hypothetical protein CDS as simple human readble tab separated values
    - `<samplename>_hypothetical.tsv`: hypothetical protein CDS amino acid sequences as FASTA
    - `<samplename>.tsv`: annotations as simple human readble TSV
    - `<samplename>.txt`: summary in TXT format

> Descriptions taken from the [Bakta documentation](https://github.com/oschwengers/bakta#output).

</details>

[Bakta](https://github.com/oschwengers/bakta) is a tool for the rapid & standardised annotation of bacterial genomes and plasmids from both isolates and MAGs. It provides dbxref-rich, sORF-including and taxon-independent annotations in machine-readable JSON & bioinformatics standard file formats for automated downstream analysis. The output is used by some of the functional screening tools.

### Protein annotation

[InterProScan](#interproscan)

#### InterProScan

<details markdown="1">
<summary>Output files</summary>

- `interproscan/`
  - `<samplename>_cleaned.faa`: clean version of the fasta files (in amino acid format) generated by one of the annotation tools (i.e. Pyrodigal, Prokka, Bakta). These contain sequences with no special characters (for eg. `*` or `-`).
  - `<samplename>_interproscan_faa.tsv`: predicted proteins and domains using the InterPro database in TSV format

</details>

[InterProScan](https://github.com/ebi-pf-team/interproscan) is designed to predict protein functions and provide possible domain and motif information of the coding regions. It utilizes the InterPro database that consists of multiple sister databases such as PANTHER, ProSite, Pfam, etc. More details can be found in the [documentation](https://interproscan-docs.readthedocs.io/en/latest/index.html).

### ARG detection tools

[AMRFinderPlus](#amrfinderplus), [hAMRonization](#hamronization), [argNorm](#argnorm)

#### AMRFinderPlus

<details markdown="1">
<summary>Output files</summary>

- `amrfinderplus/`
  - `*.tsv`: search results in tabular format

</details>

[AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder) relies on NCBI’s curated Reference Gene Database and curated collection of Hidden Markov Models. It identifies antimicrobial resistance genes, resistance-associated point mutations, and select other classes of genes using protein annotations and/or assembled nucleotide sequences.

#### hAMRonization

<details markdown="1">
<summary>Output files</summary>

- `hamronization_summarize/` one of the following:
  - `hamronization_combined_report.json`: summarised output in .json format
  - `hamronization_combined_report.tsv`: summarised output in .tsv format when the taxonomic classification is turned off (pipeline default).
  - `hamronization_complete_summary_taxonomy.tsv.gz`: summarised output in gzipped format when the taxonomic classification is turned on by `--run_taxa_classification`.
  - `hamronization_combined_report.html`: interactive output in .html format

</details>
<details markdown="1">
<summary>ARG summary table headers</summary>

| Table column                                      | Description                                                                                                                                                                         |
| ------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `input_file_name`                                 | Name of the file containing the sequence data to be analysed                                                                                                                        |
| `gene_symbol`                                     | Short name of a gene; a single word that does not contain white space characters. It is typically derived from the gene name                                                        |
| `gene_name`                                       | Name of a gene                                                                                                                                                                      |
| `reference_database_name`                         | Identifier of a biological or bioinformatics database                                                                                                                               |
| `reference_database_version`                      | Version of the database containing the reference sequences used for analysis                                                                                                        |
| `reference_accession`                             | Identifier that specifies an individual sequence record in a public sequence repository                                                                                             |
| `analysis_software_name`                          | Name of a computer package, application, method or function used for the analysis of data                                                                                           |
| `analysis_software_version`                       | Version of software used to analyze data                                                                                                                                            |
| `genetic_variation_type`                          | Class of genetic variation detected                                                                                                                                                 |
| `antimicrobial_agent` (optional)                  | A substance that kills or slows the growth of microorganisms, including bacteria, viruses, fungi and protozoans                                                                     |
| `coverage_percentage` (optional)                  | Percentage of the reference sequence covered by the sequence of interest                                                                                                            |
| `coverage_depth` (optional)                       | Average number of reads representing a given nucleotide in the reconstructed sequence                                                                                               |
| `coverage_ratio` (optional)                       | Ratio of the reference sequence covered by the sequence of interest.                                                                                                                |
| `drug_class` (optional)                           | Set of antibiotic molecules, with similar chemical structures, molecular targets, and/or modes and mechanisms of action                                                             |
| `input_gene_length` (optional)                    | Length (number of positions) of a target gene sequence submitted by a user                                                                                                          |
| `input_gene_start` (optional)                     | Position of the first nucleotide in a gene sequence being analysed (input gene sequence)                                                                                            |
| `input_gene_stop` (optional)                      | Position of the last nucleotide in a gene sequence being analysed (input gene sequence)                                                                                             |
| `input_protein_length` (optional)                 | Length (number of positions) of a protein target sequence submitted by a user                                                                                                       |
| `input_protein_start` (optional)                  | Position of the first amino acid in a protein sequence being analysed (input protein sequence)                                                                                      |
| `input_protein_stop` (optional)                   | Position of the last amino acid in a protein sequence being analysed (input protein sequence)                                                                                       |
| `input_sequence_id` (optional)                    | Identifier of molecular sequence(s) or entries from a molecular sequence database                                                                                                   |
| `nucleotide_mutation` (optional)                  | Nucleotide sequence change(s) detected in the sequence being analysed compared to a reference                                                                                       |
| `nucleotide_mutation_interpretation` (optional)   | Description of identified nucleotide mutation(s) that facilitate clinical interpretation                                                                                            |
| `predicted_phenotype` (optional)                  | Characteristic of an organism that is predicted rather than directly measured or observed                                                                                           |
| `predicted_phenotype_confidence_level` (optional) | Confidence level in a predicted phenotype                                                                                                                                           |
| `amino_acid_mutation` (optional)                  | Amino acid sequence change(s) detected in the sequence being analysed compared to a reference                                                                                       |
| `amino_acid_mutation_interpretation` (optional)   | Description of identified amino acid mutation(s) that facilitate clinical interpretation.                                                                                           |
| `reference_gene_length` (optional)                | Length (number of positions) of a gene reference sequence retrieved from a database                                                                                                 |
| `reference_gene_start` (optional)                 | Position of the first nucleotide in a reference gene sequence                                                                                                                       |
| `reference_gene_stop` (optional)                  | Position of the last nucleotide in a reference sequence                                                                                                                             |
| `reference_protein_length` (optional)             | Length (number of positions) of a protein reference sequence retrieved from a database                                                                                              |
| `reference_protein_start` (optional)              | Position of the first amino acid in a reference protein sequence                                                                                                                    |
| `reference_protein_stop` (optional)               | Position of the last amino acid in a reference protein sequence                                                                                                                     |
| `resistance_mechanism` (optional)                 | Antibiotic resistance mechanisms evolve naturally via natural selection through random mutation, but it could also be engineered by applying an evolutionary stress on a population |
| `strand_orientation` (optional)                   | Orientation of a genomic element on the double-stranded molecule                                                                                                                    |
| `sequence_identity` (optional)                    | Sequence identity is the number (%) of matches (identical characters) in positions from an alignment of two molecular sequences                                                     |

</details>

[hAMRonization](https://github.com/pha4ge/hAMRonization) summarizes the output of **AMRFinderPlus** antimicrobial resistance gene detection into a standardized unified tabular format. It supports a variety of summary options including an interactive summary.

#### argNorm

<details markdown="1">
<summary>Output files</summary>

- `normalized/`
  - `*.{tsv}`: search results in tabular format
  </details>
  <details markdown="1">
  <summary>ARG summary table headers</summary>

| Table column                 | Description                                                                      |
| ---------------------------- | -------------------------------------------------------------------------------- |
| `ARO`                        | ARO accessions of ARG                                                            |
| `confers_resistance_to`      | ARO accessions of drugs to which ARGs confer resistance to                       |
| `resistance_to_drug_classes` | ARO accessions of drugs classes to which drugs in `confers_resistance_to` belong |

</details>

[argnorm](https://github.com/BigDataBiology/argNorm) is a tool to normalize antibiotic resistance genes (ARGs) by mapping them to the antibiotic resistance ontology ([ARO](https://obofoundry.org/ontology/aro.html)) created by the CARD database. argNorm also enhances antibiotic resistance gene annotations by providing categorization of the drugs that antibiotic resistance genes confer resistance to.

argNorm takes the output of [hAMRonization](#hamronization) from [AMRFinderPlus](#amrfinderplus) and normalizes ARGs to the ARO.

### Input contig QC

<details markdown="1">
<summary>Output files</summary>

- `seqkit/`
  - `<samplename>.fasta`: FASTA file containing filtered contigs
  </details>

[SeqKit](https://bioinf.shenwei.me/seqkit/) is a cross-platform and ultrafast toolkit for FASTA/Q file manipulation used for sequence quality control and filtering.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
