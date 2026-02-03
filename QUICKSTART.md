# AMRfinderflow Quick Start Guide

## What This Pipeline Does

AMRfinderflow is a streamlined pipeline for detecting antibiotic resistance genes (ARGs) in bacterial genomes using:
- **AMRFinderPlus** for ARG detection
- **argNorm** for standardizing results to the Antibiotic Resistance Ontology

## Quick Start

### 1. Prepare Your Input

Create a samplesheet (`samplesheet.csv`):
```csv
sample,fasta
sample1,/path/to/sample1.fasta
sample2,/path/to/sample2.fasta
sample3,/path/to/sample3.fasta
```

### 2. Run the Pipeline

**Basic usage (ARG screening only):**
```bash
nextflow run barbarahelena/amrfinderflow \
  -profile docker \
  --input samplesheet.csv \
  --outdir results
```

**With taxonomic classification:**
```bash
nextflow run barbarahelena/amrfinderflow \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --run_taxa_classification true \
  --taxa_classification_mmseqs_db /path/to/mmseqs2/db
```

**With protein annotation:**
```bash
nextflow run barbarahelena/amrfinderflow \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --run_protein_annotation true
```

**Full pipeline with all features:**
```bash
nextflow run barbarahelena/amrfinderflow \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --run_taxa_classification true \
  --run_protein_annotation true \
  --annotation_tool pyrodigal
```

## Key Parameters

### Required
- `--input` - Path to samplesheet CSV file
- `--outdir` - Output directory for results
- `-profile` - Execution profile (docker, singularity, conda, etc.)

### Optional but Recommended
- `--run_taxa_classification` - Enable taxonomic classification (default: false)
- `--run_protein_annotation` - Enable InterProScan annotation (default: false)
- `--annotation_tool` - Choose annotation tool (default: pyrodigal)
  - Options: `prodigal`, `pyrodigal`, `prokka`, `bakta`

### ARG-Specific Options
- `--arg_amrfinderplus_db` - Path to custom AMRFinderPlus database
- `--arg_skip_argnorm` - Skip argNorm normalization (default: false)
- `--arg_amrfinderplus_identmin` - Minimum identity for AMRFinderPlus (default: -1)
- `--arg_amrfinderplus_coveragemin` - Minimum coverage for AMRFinderPlus (default: 0.5)

## Expected Output

```
results/
├── pipeline_info/
│   └── nf_core_funcscan_software_versions.yml
├── seqkit/
│   └── stats/
├── annotation/
│   ├── pyrodigal/  (or prodigal/prokka/bakta)
│   └── ...
├── amrfinderplus/
│   ├── sample1.tsv
│   ├── sample2.tsv
│   └── ...
├── hamronization/
│   ├── amrfinderplus/
│   └── summary/
│       └── hamronization_combined_report.tsv
└── argnorm/
    ├── sample1_normed.tsv
    ├── sample2_normed.tsv
    └── ...
```

### Key Output Files

1. **AMRFinderPlus results**: `amrfinderplus/*.tsv`
   - Raw ARG predictions from AMRFinderPlus
   
2. **argNorm normalized results**: `argnorm/*_normed.tsv`
   - ARGs mapped to ARO ontology
   
3. **HAMRONIZATION summary**: `hamronization/summary/hamronization_combined_report.tsv`
   - Standardized report combining all samples

4. **Annotation files**: `annotation/[tool]/`
   - GBK, FAA, and FNA files from genome annotation

## Profile Options

- `docker` - Use Docker containers (recommended)
- `singularity` - Use Singularity containers
- `conda` - Use Conda environments
- `test` - Run test dataset (small, quick)

## Common Use Cases

### 1. Screen for ARGs only (fastest)
```bash
nextflow run barbarahelena/amrfinderflow \
  -profile docker \
  --input samples.csv \
  --outdir results
```

### 2. ARGs with taxonomic context
```bash
nextflow run barbarahelena/amrfinderflow \
  -profile docker \
  --input samples.csv \
  --outdir results \
  --run_taxa_classification true
```

### 3. Complete analysis with annotation
```bash
nextflow run barbarahelena/amrfinderflow \
  -profile docker \
  --input samples.csv \
  --outdir results \
  --run_taxa_classification true \
  --run_protein_annotation true \
  --annotation_tool bakta \
  --annotation_bakta_db /path/to/bakta/db
```

### 4. Pre-annotated genomes
If you already have GBK files:
```csv
sample,fasta,faa,gbk
sample1,,,/path/to/sample1.gbk
sample2,,,/path/to/sample2.gbk
```

Then run:
```bash
nextflow run barbarahelena/amrfinderflow \
  -profile docker \
  --input preannotated.csv \
  --outdir results
```

## Troubleshooting

### Pipeline fails with "no hits found"
- Check that your input FASTA files are valid
- Ensure contigs are long enough (>= 1 kb recommended)

### AMRFinderPlus database issues
- The pipeline automatically downloads the latest AMRFinderPlus database
- To use a specific version: `--arg_amrfinderplus_db /path/to/db`

### Memory issues
- Reduce parallel processes in nextflow.config
- Use `-profile singularity` instead of docker if available
- Skip protein annotation if not needed

### Annotation tool selection
- **Pyrodigal** (default) - Fast, works with all downstream tools
- **Prodigal** - Classic tool, but may have GBK format issues
- **Prokka** - Comprehensive, slower, good annotations
- **Bakta** - Most comprehensive, requires database download

## Getting Help

1. Check the [nf-core/funcscan documentation](https://nf-co.re/funcscan/docs)
2. Review tool-specific documentation:
   - [AMRFinderPlus](https://github.com/ncbi/amr)
   - [argNorm](https://github.com/BigDataBiology/argNorm)
3. Check pipeline issues on GitHub

## Citation

If you use this pipeline, please cite:
- nf-core/funcscan: [10.5281/zenodo.7643099](https://doi.org/10.5281/zenodo.7643099)
- AMRFinderPlus: [10.1038/s41598-021-91456-0](https://doi.org/10.1038/s41598-021-91456-0)
- argNorm: [Check publication]
