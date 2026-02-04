# AMRfinderflow 

<p align="center">
  <img src="docs/images/amrfinderflow.png" alt="AMRfinderflow Logo" width="600">
</p>

This pipeline uses code and infrastructure developed and maintained by the nf-core￼community, reused here under the MIT license. 

## Introduction

**AMRfinderflow** is a simplified bioinformatics pipeline for detecting antibiotic resistance genes, based on the nf-core/funcscan pipeline. This pipeline focuses solely on ARG detection using AMRFinderPlus and argNorm, with supporting tools for quality control, taxonomic classification, genome annotation, and protein annotation.

## Pipeline summary

1. Quality control of input sequences with [`SeqKit`](https://bioinf.shenwei.me/seqkit/)
2. Taxonomic classification of contigs of **prokaryotic origin** with [`MMseqs2`](https://github.com/soedinglab/MMseqs2)
3. Annotation of assembled prokaryotic contigs with [`Prodigal`](https://github.com/hyattpd/Prodigal), [`Pyrodigal`](https://github.com/althonos/pyrodigal), [`Prokka`](https://github.com/tseemann/prokka), or [`Bakta`](https://github.com/oschwengers/bakta)
4. Annotation of coding sequences from 3. to obtain general protein families and domains with [`InterProScan`](https://github.com/ebi-pf-team/interproscan)
5. Screening contigs for antibiotic resistant gene-like sequences with [`AMRFinderPlus`](https://github.com/ncbi/amr) (protein mode with Pyrodigal); [`argNorm`](https://github.com/BigDataBiology/argNorm) is used to map the outputs of `AMRFinderPlus` to the [`Antibiotic Resistance Ontology`](https://www.ebi.ac.uk/ols4/ontologies/aro) for consistent ARG classification terms.
6. **Optional:** ARG abundance quantification workflow (when FASTQ files provided):
   - Extract nucleotide sequences of detected ARGs with metadata
   - Merge and deduplicate ARG sequences per sample group with [`CD-HIT-EST`](https://github.com/weizhongli/cdhit) (95% identity)
   - Create group-specific ARG catalogs
   - Map metagenomic reads to ARG catalogs with [`BWA-MEM2`](https://github.com/bwa-mem2/bwa-mem2)
   - Calculate abundance metrics (RPKM, RPK, Coverage, Prevalence)
   - Merge results per group and population-wide

**Note:** This simplified pipeline removes AMP screening, BGC screening, and MultiQC reporting functionality from the original nf-core/funcscan pipeline. ARG screening is enabled by default and uses only AMRFinderPlus (other ARG tools like DeepARG, fARGene, RGI, and ABRicate have been removed).


## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fasta
CONTROL_REP1,AEG588A1_001.fasta
CONTROL_REP2,AEG588A1_002.fasta
CONTROL_REP3,AEG588A1_003.fasta
```

Each row represents a (multi-)fasta file of assembled contig sequences.

**Optional:** For ARG abundance quantification via read mapping, you can also provide sample groups and FASTQ files:

`samplesheet.csv` (with groups):

```csv
sample,group,fasta
SAMPLE1_T0,0,sample1_t0.fasta
SAMPLE2_T0,0,sample2_t0.fasta
SAMPLE1_T1,1,sample1_t1.fasta
SAMPLE2_T1,1,sample2_t1.fasta
```

`samplesheet_fastqs.csv`:

```csv
sample,group,fastq_1,fastq_2
SAMPLE1_T0,0,sample1_t0_R1.fq.gz,sample1_t0_R2.fq.gz
SAMPLE2_T0,0,sample2_t0_R1.fq.gz,sample2_t0_R2.fq.gz
SAMPLE1_T1,1,sample1_t1_R1.fq.gz,sample1_t1_R2.fq.gz
SAMPLE2_T1,1,sample2_t1_R1.fq.gz,sample2_t1_R2.fq.gz
```

Samples with the same `group` ID will have their ARG sequences merged and deduplicated together before read mapping.

Now, you can run the pipeline using:

```bash
nextflow run barbarahelena/amrfinderflow \
   -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

**For ARG abundance quantification, add the FASTQ input:**

```bash
nextflow run barbarahelena/amrfinderflow \
   -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
   --input samplesheet.csv \
   --input_fastqs samplesheet_fastqs.csv \
   --outdir <OUTDIR>
```

> [!NOTE]
> ARG screening is **enabled by default**. The pipeline will automatically:
> - Run genome annotation (default: Pyrodigal)
> - Screen for ARGs using AMRFinderPlus
> - Normalize ARG results using argNorm
> 
> Optional features (disabled by default):
> - Taxonomic classification: add `--run_taxa_classification true`
> - Protein annotation with InterProScan: add `--run_protein_annotation true`
> - To disable argNorm: add `--arg_skip_argnorm true`

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/funcscan/usage) and the [parameter documentation](https://nf-co.re/funcscan/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/funcscan/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/funcscan/output).

## Credits

nf-core/funcscan was originally written by Jasmin Frangenberg, Anan Ibrahim, Louisa Perelo, Moritz E. Beber, James A. Fellows Yates.

We thank the following people for their extensive assistance in the development of this pipeline:

Adam Talbot, Alexandru Mizeranschi, Hugo Tavares, Júlia Mir Pedrol, Martin Klapper, Mehrdad Jaberi, Robert Syme, Rosa Herbst, Vedanth Ramji, @Microbion.

## Citations

If you use this pipeline, please credit nf-core/funcscan for your analysis, and cite it using the following doi: [10.5281/zenodo.7643099](https://doi.org/10.5281/zenodo.7643099)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
