# AMRfinderflow: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

AMRfinderflow is a pipeline for efficient and parallelised screening of long nucleotide sequences such as contigs for antimicrobial resistance genes. It can additionally identify the taxonomic origin of the sequences and provide protein domain annotations.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run barbarahelena/amrfinderflow --input samplesheet.csv --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

ARG screening is enabled by default. You can optionally enable taxonomic classification and/or protein annotation by adding the respective flag(s) to the command:

- `--run_taxa_classification` (for optional taxonomic annotations)
- `--run_protein_annotation` (for optional protein family and domain annotation)

For the taxonomic classification, MMseqs2 is currently the only tool implemented in the pipeline. Likewise, InterProScan is the only tool for protein sequence annotation.

**Example:** You want to run ARG screening with taxonomic classification and protein annotation:

```bash
nextflow run barbarahelena/amrfinderflow --input samplesheet.csv --outdir <OUTDIR> -profile docker --run_taxa_classification --run_protein_annotation
```

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing temporary files required for the run
<OUTDIR>        # Final results (location specified with --outdir)
.nextflow_log   # Log file from nextflow

# Other nextflow hidden files, eg. history of pipeline runs and old logs
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run barbarahelena/amrfinderflow -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Samplesheet input

barbarahelena/amrfinderflow takes FASTA files as input, typically contigs or whole genome sequences. To supply these to the pipeline, you will need to create a samplesheet with information about the samples you would like to analyse. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

The input samplesheet has to be a comma-separated file (`.csv`) with 2-5 columns and a header row as shown in the examples below.

**Required columns:**
- `sample`: Sample name
- `fasta`: Path to FASTA file (contigs/genome)

**Optional columns:**
- `group`: Group identifier for samples (see ARG workflow section below)
- `protein`: Pre-generated protein FASTA file (`.faa`)
- `gbk`: Pre-generated Genbank annotation file (`.gbk` or `.gbff`)

If you already have annotated contigs with peptide sequences and an annotation file in Genbank format, you can supply these to the pipeline using the optional `protein` and `gbk` columns. If these additional columns are supplied, pipeline annotation (i.e. with bakta, prodigal, pyrodigal or prokka) will be skipped and your corresponding annotation files used instead.

### Basic samplesheet (without pre-annotated data):

```csv title="samplesheet.csv"
sample,fasta
sample_1,/<path>/<to>/wastewater_metagenome_contigs_1.fasta.gz
sample_2,/<path>/<to>/wastewater_metagenome_contigs_2.fasta.gz
```

### With grouping for ARG read mapping:

```csv title="samplesheet.csv"
sample,group,fasta
sample_1,0,/<path>/<to>/wastewater_metagenome_contigs_1.fasta.gz
sample_2,0,/<path>/<to>/wastewater_metagenome_contigs_2.fasta.gz
sample_3,1,/<path>/<to>/wastewater_metagenome_contigs_3.fasta.gz
sample_4,1,/<path>/<to>/wastewater_metagenome_contigs_4.fasta.gz
```

### With pre-annotated data:

```csv title="samplesheet.csv"
sample,fasta,protein,gbk
sample_1,/<path>/<to>/wastewater_metagenome_contigs_1.fasta.gz,/<path>/<to>/wastewater_metagenome_contigs_1.faa,/<path>/<to>/wastewater_metagenome_contigs_1.fasta.gbk
sample_2,/<path>/<to>/wastewater_metagenome_contigs_2.fasta.gz,/<path>/<to>/wastewater_metagenome_contigs_2.faa,/<path>/<to>/wastewater_metagenome_contigs_2.fasta.gbk
```

| Column    | Description                                                                                                                                                                                                           |
| --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This will be used to name all output files from the pipeline. Spaces in sample names are automatically converted to underscores (`_`).                                                            |
| `group`   | **Optional.** Group identifier (integer or string) for grouping samples in ARG workflow. Samples with the same group ID will have their ARG sequences merged and deduplicated together before read mapping. If not provided, each sample is treated as its own group. |
| `fasta`   | Path or URL to a gzipped or uncompressed FASTA file. Accepted file suffixes are: `.fasta`, `.fna`, or `.fa`, or any of these with `.gz`, e.g. `.fa.gz`.                                                               |
| `protein` | **Optional.** Path to a pre-generated amino acid FASTA file (`.faa`) containing protein annotations of `fasta`, optionally gzipped. Required to be supplied if `gbk` also given.                                           |
| `gbk`     | **Optional.** Path to a pre-generated annotation file in Genbank format (`.gbk`, or `.gbff`) format containing annotations information of `fasta`, optionally gzipped. Required to be supplied if `protein` is also given. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### FASTQ input for ARG read mapping

If you want to perform ARG abundance quantification by mapping metagenomic reads to the detected ARG catalog, you can provide an additional FASTQ samplesheet via the `--input_fastqs` parameter:

```bash
--input_fastqs '[path to fastq samplesheet file]'
```

The FASTQ samplesheet should be a comma-separated file (`.csv`) with 4 columns (`sample`, `group`, `fastq_1`, `fastq_2`):

```csv title="samplesheet_fastqs.csv"
sample,group,fastq_1,fastq_2
sample_1,0,/<path>/<to>/sample_1_R1.fq.gz,/<path>/<to>/sample_1_R2.fq.gz
sample_2,0,/<path>/<to>/sample_2_R1.fq.gz,/<path>/<to>/sample_2_R2.fq.gz
sample_3,1,/<path>/<to>/sample_3_R1.fq.gz,/<path>/<to>/sample_3_R2.fq.gz
sample_4,1,/<path>/<to>/sample_4_R1.fq.gz,/<path>/<to>/sample_4_R2.fq.gz
```

| Column    | Description                                                                                                                                                                                                           |
| --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Sample name. Should match the `sample` names in the main input samplesheet.                                                                                                                                           |
| `group`   | Group identifier. Should match the `group` values in the main input samplesheet.                                                                                                                                      |
| `fastq_1` | Path to forward read FASTQ file (gzipped or uncompressed).                                                                                                                                                            |
| `fastq_2` | Path to reverse read FASTQ file (gzipped or uncompressed).                                                                                                                                                            |

The pipeline will:
1. Detect ARGs in contigs using AMRFinderPlus (protein mode with Pyrodigal)
2. Extract ARG nucleotide sequences with metadata (gene symbol, class, subclass)
3. Merge and deduplicate ARG sequences per group (95% identity with CD-HIT-EST)
4. Create BWA-MEM2 index for each group's ARG catalog
5. Map reads from each sample to their group's ARG catalog
6. Calculate abundance metrics (RPKM, RPK, Coverage, Prevalence) for each ARG
7. Merge results per group and across all samples

:::danger
We highly recommend performing quality control on input contigs before running the pipeline. You may not receive results for some tools if none of the contigs in a FASTA file reach certain thresholds. Check parameter documentation for relevant minimum contig parameters.
:::

## Notes on screening tools, taxonomic and functional classifications

The implementation of some tools in the pipeline may have some particular behaviours that you should be aware of before you run the pipeline.

### MMseqs2

MMseqs2 is currently the only taxonomic classification tool used in the pipeline to assign a taxonomic lineage to the input contigs. The database used to assign the taxonomic lineage can either be:

- A custom based database created by the user using `mmseqs createdb` externally and beforehand. If this flag is assigned, this database takes precedence over the default database in `--mmseqs_db_id`.

  ```bash
  --taxa_classification_mmseqs_db '<path>/<to>/<mmsesqs_custom_database>/<directory>'
  ```

  The contents of the directory should have files such as `<dbname>.version` and `<dbname>.taxonomy` in the top level.

- An MMseqs2 ready database. These databases were compiled by the developers of MMseqs2 and can be called using their labels. All available options can be found [here](https://github.com/soedinglab/MMseqs2/wiki#downloading-databases). Only use those databases that have taxonomy files available (i.e. Taxonomy column shows "yes"). By default MMseqs2 in the pipeline uses '[Kalamari](https://github.com/lskatz/Kalamari)', and runs an amino acid-based alignment. However, if the user requires a more comprehensive taxonomic classification, we recommend the use of [GTDB](https://gtdb.ecogenomic.org/), but for that please remember to increase the memory, CPU threads and time required for the process `MMSEQS_TAXONOMY`.

  ```bash
  --taxa_classification_mmseqs_db_id 'Kalamari'
  ```

### InterProScan

[InterProScan](https://github.com/ebi-pf-team/interproscan) is currently the only protein annotation tool in this pipeline that gives a snapshot of the protein families and domains for each coding region.

The protein annotation workflow is activated with the flag `--run_protein_annotation`.
InterProScan is used as the only protein annotation tool at the moment and the [InterPro database](http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.72-103.0) version 5.72-103.0 is downloaded and prepared to screen the input sequences against it.

Since the database download is huge (5.5GB) and might take quite some time, you can skip the automatic database download (see section [Databases and reference files](usage/#interproscan-1) for details).

:::info
By default, the databases used by InterProScan is set as `PANTHER,ProSiteProfiles,ProSitePatterns,Pfam`. An addition of other application to the list does not guarantee that the results will be integrated correctly within `AMPcombi`.
:::

## Databases and reference files

Various tools of AMRFinderFlow use databases and reference files to operate.

nf-core/funcscan offers the functionality to auto-download databases for you, and as these databases can be very large, we suggest to store these files in a central place from where you can reuse them across pipeline runs.

If your infrastructure has internet access (particularly on compute nodes), we **highly recommend** allowing the pipeline to download these databases for you on a first run, saving these to your results directory with `--save_db`, then moving these to a different location (in case you wish to delete the results directory of this first run). An exception to this is HMM files where no auto-downloading functionality is possible.

:::warning
We generally do not recommend downloading the databases yourself, as this can often be non-trivial to do!
:::

As a reference, we will describe below where and how you can obtain databases and reference files used for tools included in the pipeline.

### Bakta

AMRFinderFlow offers multiple tools for annotating input sequences. Bakta is a new tool touted as a bacteria-only successor to the well-established Prokka.

To supply the preferred Bakta database (and not have the pipeline download it for every new run), use the flag `--annotation_bakta_db`.
The full or light Bakta database must be downloaded from the Bakta Zenodo archive.

You can do this by installing via conda and using the dedicated command

```bash
conda create -n bakta -c bioconda bakta
conda activate bakta

bakta_db download --output <LOCATION_TO_STORE> --type <full|light>
```

Alternatively, you can manually download the files via the links which can be found on the [Bakta GitHub repository](https://github.com/oschwengers/bakta#database-download).

Once downloaded this must be untarred:

```bash
tar xvzf db.tar.gz
```

And then passed to the pipeline with:

```bash
--annotation_bakta_db /<path>/<to>/<db>/
```

The contents of the directory should have files such as `*.dmnd` in the top level.

:::info
The flag `--save_db` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

### AMRFinderPlus

AMRFinderPlus relies on NCBI's curated Reference Gene Database and curated collection of Hidden Markov Models.

AMRfinderflow will download this database for you, unless the path to a local version is given with:

```bash
--arg_amrfinderplus_db '/<path>/<to>/<amrfinderplus_db>/latest'
```

You must give the `latest` directory to the pipeline, and the contents of the directory should include files such as `*.nbd`, `*.nhr`, `versions.txt` etc. in the top level.

To obtain a local version of the database:

1. Install AMRFinderPlus from [bioconda](https://bioconda.github.io/recipes/ncbi-amrfinderplus/README.html?highlight=amrfinderplus).
   To ensure database compatibility, please use the same version as is used in your amrfinderflow release (check version in file `<installation>/<path>/amrfinderflow/modules/nf-core/amrfinderplus/run/environment.yml`).

```bash
conda create -n amrfinderplus -c bioconda ncbi-amrfinderplus=3.12.8
conda activate amrfinderplus
```

2. Run `amrfinder --update`, which will download the latest version of the AMRFinderPlus database to the default location (location of the AMRFinderPlus binaries/data).
   It creates a directory in the format YYYY-MM-DD.version (e.g., `<installation>/<path>/data/2024-01-31.1/`).

<details markdown="1">
<summary>AMR related files in the database folder</summary>

```tree
<YYYY-MM-DD.v>/
├── AMR_CDS.*
├── AMR_DNA-Campylobacter.*
├── AMR_DNA-Clostridioides_difficile.*
├── AMR_DNA-Enterococcus_faecalis.*
├── AMR_DNA-Enterococcus_faecium.*
├── AMR_DNA-Escherichia.*
├── AMR_DNA-Neisseria.*
├── AMR_DNA-Salmonella.*
├── AMR_DNA-Staphylococcus_aureus.*
├── AMR_DNA-Streptococcus_pneumoniae.*
├── AMR.LIB.*
├── AMRProt.*
├── changes.txt
├── database_format_version.txt
├── fam.tab
├── taxgroup.tab
└── version.txt
```

</details>

:::info
The flag `--save_db` saves the pipeline-downloaded databases in your results directory. You can then move these to a central cache directory of your choice for re-use in the future.
:::

### MMSeqs2

To download MMSeqs2 databases for taxonomic classification, you can install `mmseqs` via conda:

```bash
conda create -n mmseqs2 -c bioconda mmseqs2
conda activate mmseqs2
```

Then to download the database of your choice

```bash
mmseqs databases <DATABASE_NAME> <LOCATION_TO_STORE> tmp/
```

:::info
You may want to specify a different location for `tmp/`, we just borrowed here from the official `mmseqs` [documentation](https://github.com/soedinglab/mmseqs2/wiki#downloading-databases).
:::

### InterProScan

[InterProScan](https://github.com/ebi-pf-team/interproscan) is used to provide more information about the proteins annotated on the contigs. By default, turning on this subworkflow with `--run_protein_annotation` will download and unzip the [InterPro database](http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.72-103.0/) version 5.72-103.0. The database can be saved in the output directory `<output_directory>/databases/interproscan/` if the `--save_db` is turned on.

:::note
The huge database download (5.5GB) can take up to 4 hours depending on the bandwidth.
:::

A local version of the database can be supplied to the pipeline by passing the InterProScan database directory to `--protein_annotation_interproscan_db <path/to/downloaded-untarred-interproscan_db-dir/>`. The directory can be created by running (e.g. for database version 5.72-103.0):

```
curl -L https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.72-103.0/interproscan-5.72-103.0-64-bit.tar.gz -o interproscan_db/interproscan-5.72-103.0-64-bit.tar.gz
tar -xzf interproscan_db/interproscan-5.72-103.0-64-bit.tar.gz -C interproscan_db/

```

The contents of the database directory should include the directory `data` in the top level with a couple of subdirectories:

```
interproscan_db/
    └── data/
    ├── antifam
    ├── cdd
    ├── funfam
    ├── gene3d
    ├── hamap
    ├── ncbifam
    ├── panther
    | └── [18.0]
    ├── pfam
    | └── [36.0]
    ├── phobius
    ├── pirsf
    ├── pirsr
    ├── prints
    ├── prosite
    | └── [2023_05]
    ├── sfld
    ├── smart
    ├── superfamily
    └── tmhmm
```

## Updating the pipeline

When you run the below command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull barbarahelena/amrfinderflow
```

## Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [barbarahelena/amrfinderflow releases page](https://github.com/barbarahelena/amrfinderflow/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
