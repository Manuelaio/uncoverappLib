This repository is home of R packages uncoverappLib launching a *unCOVERApp*, a web application for clinical assessment and annotation of coverage gaps in target genes. 
Read more about unCOVERApp [here](https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1)

[![DOI](https://zenodo.org/badge/254597958.svg)](https://zenodo.org/badge/latestdoi/254597958)



# Table of content

* [Prerequisites](#Prerequisites)
* [Installation](#Installation)
* [Input_file_preparation](#Input_file_preparation)
* [Usage](#Usage)


## Prerequisites


This app requires following dependencies:
- samtools v.1.9
- R v.3.5.1 or RStudio
- annotation files that can be downloaded on Zenodo at the following link https://zenodo.org/record/3747448#.XpBmnVMzbOR and stored in a user folder. 

## Installation


The development version can be installed directly from Github:

``` {r}
install.packages("devtools")
devtools::install_github("Manuelaio/uncoverappLib")
``` 

To run locally, users need to install this library and run. It is mandatory to expliciatate absolute path of prevusly downloaded annotation files before run the app. 

``` {r}
library(uncoverappLib)
file.name='../path/sorted.bed.gz'
tbi='.../path/sorted.bed.gz.tbi'
run.uncoverapp()

``` 
The rapid spread of NGS-technology has been coupled since its beginning with development of bioinformatic tools for data analysis and interpretation. However, despite increasing accuracy of available approaches, the need to assess sequencing quality of the analysis targets at the base-pair resolution poses growing challenges especially in the clinical settings.  In diagnostics indeed, meticulous investigation of every single target base is often required to exclude that pathogenic events across the gene of interest may be missed due to uneven sequence coverage.


unCOVERApp is an interactive web-application developmented for graphical inspection of sequence coverage within gene regions.


unCOVERApp highlights low coverage genomic positions, according to the coverage threshold specified by the user, providing dbNSFP-based annotations for clinical assessment of low coverage regions. It implements basic statistical tools such as binomial probability calculation that genomic positions are adequately covered, and [maximum credible allele frequency](http://cardiodb.org/allelefrequencyapp/). 


# Input_file_preparation

To associate low-coverage sites with functional and clinical annotations, unCOVERApp uses dbNSFP version 4.0 stored in two file:


* `sorted.bed.gz`: a genomically-sorted, TABIX-indexed, BGZipped BED file containing selected columns from dbNSFP version  v4.0. Only following columns from chromosome-specific dbNSFP files were be selected (in that specific order): $1,$9,$9,$3,$4,$7,$13,$16,$2,$2,$39,$48,$78,$104,$229,$365,$370,$371 and merged into a single file. The resulting file (i.e. db_all.bed) were be then sorted, BGZipped as below. The output of these commands were be redirected to file "sorted.bed.gz" and TABIX-indexed.  This file is available on Zenodo [here](https://zenodo.org/record/3747448#.XpBmnVMzbOR)


* `sorted.bed.gz.tbi`: TABIX-indexed file. This file is available on Zenodo [here](https://zenodo.org/record/3747448#.XpBmnVMzbOR)


In the home page, unCOVERApp serves following static files in order to prepare input file to load on application.  

* `preprocessing.R ` 

* `uncoverappLib.config `

* `make_bed.sh `: bash script to create the input file for unCOVERApp. 

Users must download and store three scripts in a folder. 
Compile a `configuration file ` specifying absolute path of: unCOVERApp script folder, txt file containing HGNC gene name(s)(one per row), and path to file with ".list" extension containing absolute paths to BAM files (one per row) and folder output location. Compile genome reference (hg19 or hg38) and chromosome notation BAM (number refers to 1, 2, ???X,.M notation BAM, chr refers to chr1, chr2,??? chrX, chrM notation BAM). The resulting unCOVERApp input file is a BED file (tab-separated) containing depth of coverage (DoC) for each genomic position (one per row) of target genes for as many samples as many BAM files are listed in the ".list" file. Finally you run following bash command:

```sh
bash www/make_bed.sh www/uncoverappLib.config

#substitute www with your folder name

```

[![Build Status](https://travis-ci.com/Manuelaio/test_app.svg?token=25AMAYuQwZENC1xVJVSe&branch=master)](https://travis-ci.com/Manuelaio/test_app)

The log file is a trouble shooter, so please revise when any problem happens. In order to have more information about dependeces on [github](https://github.com/Manuelaio/test_dependence)

[![Build Status](https://travis-ci.com/Manuelaio/test_dependence.svg?branch=master)](https://travis-ci.com/Manuelaio/test_dependence)

Bash script creates a new directory named with current date in users-defined location, in which input file is stored called *multisample.bed.gz file*.


# Usage

Load library and set up R environment with annotation file as following. 
The way to launchs unCOVERApp is *run.uncoverapp()* function. 

``` {r}
library(uncoverappLib)
file.name='../path/sorted.bed.gz'
tbi='.../path/sorted.bed.gz.tbi'
run.uncoverapp()

``` 

For more details on working with uncoverapp see Vignette or [Documentation.pdf](https://github.com/Manuelaio/unCOVERApp/blob/master/Documentation.pdf) on Github. 


