# Methylation analysis pipeline and scripts for the MethylNads project

### Rose Driscoll, Cynthia O'Rourke, Josh Faber-Hammond, Pete Hurd, and Suzy Renn

## Background

This repository contains the data processing pipelines and a few custom scripts associated with the MethylNads project [manuscript title to be added]

## Data

The raw data for this project is available from GEO [accession information to be added]

A metadata file (`metadata.txt`) with information about each of the samples is provided in this repository.

The analysis pipelines require *in silico* bisulfite converted reference amplicon sequences; these are provided in the `reference_sequences` directory [to be added]

## Scripts

The data was analyzed in two different ways, so two pipelines are supplied in this repository.

### Percent methylation analysis pipeline

The `percent-methylation` directory contains scripts relevant to the percent methylation analysis pipeline. This pipeline calculates percent methylation for each sample for each CpG, as well as creating associated files that may be used to check the completeness of the bisulfite conversion reaction and the sequencing error rate.

Bash scripts for this pipeline (to be run in this order):

- `methylation_pipeline_bowtie_20180605_RD.txt`
- `methylation_pipeline_samtoolspythonandcat_20180605_RD.txt`

Custom python scripts required by the pipeline:

- `basecount_20180518_RD.py` (see note about this script under "Authorship and acknowledgements")
- `extract_CpGs_20180521_RD.py`
- `conversion_completeness_20180521_RD.py`
- `A_G_error_20180521_RD.py`

Header files, also required by the pipeline:

- `head_CpGs_20180521_RD.txt`
- `head_conversion_20180521_RD.txt`
- `head_error_20180521_RD.txt`

### Epiallele analysis pipeline

The `epialleles` directory contains scripts relevant to the epiallele analysis pipeline. This pipeline calculates the frequency of each epiallele (unique combination of methylated and unmethylated CpGs) for each sample for each amplicon.

Bash scripts for this pipeline (to be run in this order):

- `flash_and_bowtie2_20180606_RD.txt`
- `AI_usearch_and_muscle_20180606_RD.txt` |
- `AP_usearch_and_muscle_20180606_RD.txt` |
- `B2_usearch_and_muscle_20180606_RD.txt` | (these four may be run in any order, or may be run concurrently)
- `BI_usearch_and_muscle_20180606_RD.txt` | 
- `end_usearch_and_muscle_20180606_RD.txt`

Custom python scripts required by the pipeline:

- `find_CpGs_AI.py`
- `find_CpGs_AP.py`
- `find_CpGs_AP6.py`
- `find_CpGs_B2.py`
- `find_CpGs_BI.py`

This pipeline creates its own header, so a header file is not supplied.

## Dependencies

The following programs are required to run the analysis pipelines:

- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [SAMtools](http://www.htslib.org/)
- [Python 2.7](https://www.python.org/)
- [USEARCH](https://www.drive5.com/usearch/)
- [MUSCLE](https://www.drive5.com/muscle/)
- [FASTX toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
- [Clustal Omega](http://www.clustal.org/omega/)

## Authorship and acknowledgements

`basecount_20180518_RD.py` was originally based on a script authored by Damian Kao which is publicly available [here](http://blog.nextgenetics.net/?e=56). All other scripts in this repository are the original work of Rose Driscoll and Josh Faber-Hammond. 