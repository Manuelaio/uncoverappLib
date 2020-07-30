### An interactive graphical web application for clinical assessment of sequence coverage at base-pair level

This is a web application for clinical assessment of sequence coverage. 
unCOVERApp allows:


- to display interactive plots showing sequence gene coverage down to base-pair
resolution and functional/ clinical annotations of sequence 
positions within coverage gaps (` Coverage Analysis` page).


- to calculate the [maximum credible population allele frequency](http://cardiodb.org/allelefrequencyapp/) (AF) to be applied as AF 
filtering threshold tailored to the model of the disease-of-interest 
instead of a general AF cut-off (e.g. 1 % or 0.1 %) 
(` Calculate AF by allele frequency app ` page).



- to calculate the 95 % probability of the binomial distribution to observe at 
least N variant-supporting reads (N is the number of successes) based on a 
user-defined allele fraction that is expected for the variant 
(which is the probability of success). Especially useful to obtain the 
range of variant-supporting reads that is most likely to occur at a given 
a depth of coverage (DoC)(which is the number of trials) for somatic variants
with low allele fraction
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
included in extdata of uncoverappLib packages

Below is an example of genes list. 


```{r}
DNAJC8
GNB1
PEX10
RPL22
```

- a text file, with .list extension, containing absolute paths to BAM files
(one per row) and to be uploaded to ` Load bam file(s) list ` box.
In the output file, sample 1,2,3.., correspond
to the sample in the bam file bam.list file listed in row 1,2,3, …. 
Below is an example of bam.list in a Unix-like OS including macOS. 

```{r}
/home/user/bam/sample1.bam
/home/user/bam/sample2.bam
/home/user/bam/smaple3.bam
```

While all inputs are loading, a progress bar appears during processing phase. 
unCOVERApp input file generation fails if incorrect gene names are specified. 
An unrecognized gene name(s) table is displayed if such a case occurs.

Below is a snippet of bed file output of **Preprocessing** 
unCOVERApp. 

```{r}

chr15   89859516        89859516        68      A:68
chr15   89859517        89859517        70      T:70
chr15   89859518        89859518        73      A:2;G:71
chr15   89859519        89859519        73      A:73
chr15   89859520        89859520        74      C:74
chr15   89859521        89859521        75      C:1;T:74

```


The preprocessing time depends on the size of the BAM file(s) and on the number 
of genes to investigate. In general, if many (e.g. > 50) genes are to be analyzed, 
we would recommend using `buildInput` function and run it in R console 
before launching  the app.

Alternatively, other tools do a similar job and can be used to generate the 
unCOVERApp input file ( for instance:
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


- ` coverage threshold ` : required coverage threshold  

- ` Sample  ` : sample(s) to analyze according to help text indications on the 
    App page

- ` Transcript number ` : transcript number as in first column 
    of ` Exon genomic coordinate positions from UCSC ` output App table.

- ` exon number ` and push ` Make exon ` : to zoom in a specific exon


Other input sections, as ` Chromosome `, ` Transcript ID `, 
` START genomic position `, ` END genomic position ` and ` Region coordinate `,
are dynamically filled. 


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
`Annotation on low-coverage positions ` . By clicking on the `download` button, 
users can save the table as spreadsheet format with certain cells colored 
according to pre-specified thresholds for AF, 
CADD, MAP-CAP, SIFT, Polyphen2, ClinVar, OMIM ID, HGVSp and HGVSc, ...).


- In **Calculate maximum credible allele frequency** page, users can set 
allele frequency cut-offs based on specific assumptions about the genetic 
architecture of the disease. Users may click on the `download` button and 
save the resulting table as spreadsheet format. 

- The **Binomial distribution** page returns the 95 % binomial probability 
distribution of the variant supporting reads on the input genomic position 
(`START genomic position` and `END genomic position`).
Users should define  the expected `allele fraction`
(the expected fraction of variant reads, probability of success) 
and `Variant reads ` (the minimum number of variant reads required by the user to 
support variant calling, number of successes).





