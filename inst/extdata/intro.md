### An interactive graphical web application for clinical assessment of sequence coverage at base-pair level

This is a web application for clinical assessment of sequence coverage. 
unCOVERApp allows:


- to display interactive plots showing sequence gene coverage down to base-pair
resolution and annotations of functional/clinical sequence 
positions within coverage gaps (` Coverage Analysis` page).


- to calculate maximum credible allele frequency 
(http://cardiodb.org/allelefrequencyapp/) for  allele frequency
threshold  setting to be used in variant filtering
 (` Calculate AF by allele frequency app ` page).



- to calculate the binomial probability that coverage, over a position of 
interest,  detect a variant that the user specifies 
the expected allelic fraction and number of variant reads - applicable 
especially to somatic variants
 (` Binomial distribution`page). 

### Documentation 

All unCOVERApp functionalities are based on the availability of a bed file 
containing tab-separated specification of genomic coordinates (chromosome, 
start position, end position),  coverage value and reference: alternate 
allele counts read for each position.
In the first page **Preprocessing**, users prepare the bed file providing 
input files containing a list of genes and a list of bam files, respectively: 



- a text file, with .txt extension, containing HGNC official gene name(s) one per 
row and to be uploaded to ` Load a gene(s) file ` box. An example file is
included in extdata of uncoverappLib packages, users can also download the
file `mygene.txt` also
[from](https://github.com/Manuelaio/unCOVERApp/blob/master/script/mygene.txt) 
github. 
Below is an example of genes list. 


```{r}
DNAJC8
GNB1
PEX10
RPL22
```

-  a text file, with .list extension, containing absolute paths to BAM files
(one per row) and to be uploaded to ` Load bam file(s) list ` box.
In the output file, sample 1,2,3.., correspond
to the sample in the bam file bam.list file listed in row 1,2,3, … . 
Below is an example of bam.list in a Unix-like OS including macOS. 

```{r}
/home/user/bam/sample1.bam
/home/user/bam/sample2.bam
/home/user/bam/smaple3.bam
```


While all inputs are loading, a progress bar appears during processing phase. 
Output generation will fail if there are incorrect HGNC official gene names in list.
An unrecognized gene name(s) table  will be displayed. 

Below is a snippet of bed file output of **Preprocessing** 
unCOVERApp. Users could find the example file [here](https://github.com/Manuelaio/unCOVERApp/blob/master/script/POLG.example.bed)

```{r}

chr15   89859516        89859516        68      A:68
chr15   89859517        89859517        70      T:70
chr15   89859518        89859518        73      A:2;G:71
chr15   89859519        89859519        73      A:73
chr15   89859520        89859520        74      C:74
chr15   89859521        89859521        75      C:1;T:74

```


The input BAM file(s) could be also processed on a local machine with the 
Rscript available on unCOVERAPP **Preprocessing** page
[here](https://github.com/Manuelaio/unCOVERApp/blob/master/script/Rpreprocessing.R) 
or with other tools (for instance: 
[bedtools](https://bedtools.readthedocs.io/en/latest/#), 
[samtools](http://www.htslib.org/doc/samtools-depth.html), 
[gatk](https://gatk.broadinstitute.org/hc/en-us)). 
In this case, users can load the file directly on
**Coverage Analysis**  page in `Select input file` box. 

Once bed file is ready, users can move to **Coverage Analysis** page for own 
analysis and push `load prepared input file` button.

To assess sequence coverage the following **input parameters** must be 
specified in the sidebar of the **Coverage Analysis** section:


- ` Reference Genome` : reference genome (hg19 or hg38) 

- ` Gene name ` and push ` Apply ` button:  HGNC official gene name 

-  ` Chromosome ` : chromosome number

- ` coverage threshold ` : required coverage threshold  

- ` Sample  ` : sample(s) to analyze according to help text indications on the 
App page

- ` Transcript number ` : transcript number as in first column 
of ` Exon genomic coordinate positions from UCSC ` output App table.

- ` exon number ` and push ` Make exon ` : to zoom in a specific exon


Other input sections, as ` Transcript ID `, ` START genomic position `, 
` END genomic position ` and ` Region coordinate `, are dynamically filled. 


unCOVERApp generates the following **outputs** : 


- unfiltered BED file in` bed file ` and the corresponding filtered dataset 
in ` Low coverage positions ` 

- UCSC gene coordinates in ` UCSC gene ` table

- exon coordinates and transcript number in ` UCSC exon ` table

- sequence gene coverage plot in ` Gene coverage `. The plot displays 
chromosome ideogram, genomic location and gene annotations from *Ensembl*  
and transcripts annotation from *UCSC*.
Related table shows the number of 
uncovered position in each exon given a user-defined transcript  number 

- plot of a specific exon, selected by sidebar  ` exon number ` , 
in ` Exon coverage `. Related table shows the number of Low coverage positions
annotate in ClinVar and with a high impact annotation. 

- dbNSFP annotation of low coverage positions can be found in  
`Annotation on low-coverage positions ` . Users may click on download button 
and save the table as spreadsheet format with certain cells colored  following 
thresholds for clinically-relevant variant parameters (gnomAD allele frequency,
CADD, MAP_CAP, SIFT, Polyphen2, ClinVar, OMIM ID, HGVSp and HGVSc, functional
impact of a variant score).



In **Calculate maximum credible allele frequency** page, users can set 
allele frequency cut-offs based on assumptions about the genetic architecture 
of the disease, if not calculated, variant allele frequency 5% will be used 
instead for filtering


The last page **Binomial distribution** returns the 95% binomial probability 
distribution of the variant read number on the input genomic position 
(` START genomic position` and ` END genomic position `). 

Users must specify as input the `allele fraction` (the expected fraction of
variant reads, probability of success) and `Variant reads ` (the minimum 
number of variant reads required by the user to support variant calling,
number of successes). 



