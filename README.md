# Bisulfite amplicon sequencing methylation analysis pipeline and scripts

### Rose Driscoll, Cynthia O'Rourke, Josh Faber-Hammond, Pete Hurd, and Suzy Renn

Corresponding author: <renns@reed.edu>

## Background

This repository contains the data processing pipelines and custom scripts associated with the paper "Epigenetic regulation of gonadal and brain aromatase expression in a cichlid fish with environmental sex determination" (Driscoll *et al* 2020).

## Data

The raw data for this project is available from GEO (#GSE135681)

The analysis pipelines require *in silico* bisulfite converted reference amplicon sequences; these are provided in the `reference_sequences` directory.

Processed data files for the samples used in the statistical analysis for Driscoll *et al* 2020 are provided in the `statistical-analysis` directory (described below).

## Data processing scripts

The sequence data was analyzed in two different ways; both pipelines are supplied in this repository and described below.

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
- `AI_usearch_and_muscle_20180606_RD.txt` \*
- `AP_usearch_and_muscle_20180606_RD.txt` \*
- `B2_usearch_and_muscle_20180606_RD.txt` \* 
- `BI_usearch_and_muscle_20180606_RD.txt` \*
- `end_usearch_and_muscle_20180606_RD.txt`

\* These four scripts may be run in any order, or may be run concurrently.

Custom python scripts required by the pipeline:

- `find_CpGs_AI.py`
- `find_CpGs_AP.py`
- `find_CpGs_AP6.py`
- `find_CpGs_B2.py`
- `find_CpGs_BI.py`

This pipeline creates its own header, so a header file is not supplied.

### Dependencies

The following programs are required to run the analysis pipelines:

- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.2.9
- [SAMtools](http://www.htslib.org/) v1.9
- [Python 2.7](https://www.python.org/)
- [USEARCH](https://www.drive5.com/usearch/) v5.2.236
- [MUSCLE](https://www.drive5.com/muscle/) v3.8.31
- [FASTX toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
- [Clustal Omega](http://www.clustal.org/omega/)

### Authorship and acknowledgements

`basecount_20180518_RD.py` was originally based on a script authored by Damian Kao which is publicly available [here](http://blog.nextgenetics.net/?e=56). All other data processing scripts are the original work of Rose Driscoll and Josh Faber-Hammond. 

## Statistical analysis

The `statistical-analysis` directory contains scripts used to perform statistical analyses and visualize (i.e., make figures), plus corresponding processed data files.

### Scripts

The following statistical analysis scripts are supplied in this directory:

- `20190802_epiallele_Driscoll_2019.Rmd`
- `expression.r`
- `pct-methyl.r`

The following visualization (figure-making) scripts are supplied in this directory:

- `expression-graphs.r`
- `pct-methyl-graphs.r`

### Processed data files

The following processed data files are supplied in this directory:

- `driscoll_2019_epiallele_data_counts.csv`
- `expression_qPCR.csv`
- `pct-methylation.csv`

In addition, a metadata file is supplied: `driscoll_2019_metadata.csv`

### Dependencies

Statistical analysis and visualization for the manuscript was performed with R version 3.6.0. The following R packages are required to run the statistical analysis and visualization scripts:

- [dplyr](https://CRAN.R-project.org/package=dplyr)
- [reshape](https://CRAN.R-project.org/package=reshape)
- [car](https://CRAN.R-project.org/package=car)
- [lme4](https://CRAN.R-project.org/package=lme4)
- [lmerTest](https://CRAN.R-project.org/package=lmerTest)
- [emmeans](https://CRAN.R-project.org/package=emmeans)

### Authorship

All statistical analysis and visualization scripts are the original work of Suzy Renn and Pete Hurd. 

