### An interactive graphical web application for clinical assessment of sequence coverage at base-pair level

This is a web application for graphical inspection of sequence coverage within gene regions. 
unCOVERApp allows:

- interactive plots displaying sequence gene coverage down to base-pair resolution and functional and clinical downloadable annotation of base-pair  positions within user-defined coverage gaps (` Coverage Analysis` page).


- calculation and definition of maximum credible allele frequency (http://cardiodb.org/allelefrequencyapp/) to be used as allele frequency filtering threshold (` Calculate AF by allele frequency app ` page).

- calculation of binomial probability, for somatic variants, that position coverage is adequate to detect variant allele provided the expected allelic fraction and the required number of variant reads  (` Binomial distribution` page). 

### Documentation 

The input BAM file(s) should be processed on a local machine with the bash script available in GitHub [here](https://github.com/Manuelaio/unCOVERApp/blob/master/www/make_bed.sh) and also at the end of that document. Downloading unCOVERApp folder from GitHub, users found this preprocessing scritp in ` www` folder. 


The bash script requires R v.3.5.1, samtools v.1.9, the companion Rscript ` preprocessing.R ` available in GitHub [here](https://github.com/Manuelaio/unCOVERApp/blob/master/www/preprocessing.R) and a easy  [configuration file](https://github.com/Manuelaio/unCOVERApp/blob/master/www/uncoverappLib.config) to function. Users instruction is available in GitHub on [README.md](https://github.com/Manuelaio/unCOVERApp/blob/master/README.md)



BAM file processing will generate a BED-gzipped file (.bed.gz) containing per-sample base coverage counts.


Once uploaded to unCOVERApp, the BED file content can be visualized in table ` bed file `.


To assess sequence coverage the following **input parameters** must be specified in the sidebar of the **Coverage Analysis** section:


- ` Reference Genome` : reference genome (hg19 or hg38) 

- ` Gene name ` and push ` Apply ` button:  HGNC official gene name 

-  ` Chromosome ` : chromosome number

- ` coverage threshold ` : required coverage threshold  

- ` Sample  ` : sample(s) to analyze according to help text indications on the App page

- ` Transcript number ` : transcript number as in first column of ` Exon genomic coordinate positions from UCSC ` output App table.

- ` exon number ` and push ` Make exon ` : to zoom in a specific exon


Other input sections, as ` Transcritp ID `, ` START genomic position `, ` END genomic position ` and ` Region coordinate `, are dynamically filled. 


unCOVERApp generates the following **outputs** : 


- unfiltered BED file in` bed file ` and the corresponding filtered dataset in ` Low coverage positions ` 

- sequence gene coverage plot in ` Gene coverage `. The plot displays chromosome ideogram, genomic location and gene annotations from *Ensembl*  and transcripts annotation from *UCSC*. Related table shows the number of uncovered position in each exons given a user-defined transcript  number 

- plot of a specific exon, selected by sidebar  ` exon number ` , in ` Exon coverage `. Related table shows the number of Low coverage positions present in ClinVar and with a high impact annotation. 

- dbNSFP annotation of low coverage positions can be found in  `Annotation on low-coverage positions ` . Users may click on download button and save the table as spreadsheet format with certain cells colored  following thresholds for clinically-relevant variant parameters (gnomAD allele frequency, CADD, MAP_CAP, SIFT, Polyphen2, Clinvar, OMIM ID, HGVSp and HGVSc, functional impact of a variant score).

- UCSC gene coordinates in ` UCSC gene ` table

- exon coordinates and transcript number in ` UCSC exon ` table


In **Calculate maximum credible allele frequency** page, users can set allele frequency cut-offs based on assumptions about the genetic architecture of the disease, if not calculated, variant allele frequency 5% will be used instead for filtering


The last page **Binomial distribution** returns the 95% binomial probability distribution of the variant read number on the input genomic position (` START genomic position` and ` END genomic position `). 

Users must specify as input the `allele fraction` (the expected fraction of variant reads, probability of success) and `Variant reads ` (the minimum number of variant reads required by the user to support variant calling, number of successes). 



