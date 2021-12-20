### Preparation input file 
In this section users can prepare bed file required for clinical 
assessment of low depth of coverage genomic regions. 
The resulting unCOVERApp input file is a BED file (tab-separated) containing a 
minimum of five columns, chromosome, start and end position, depth of coverage 
(DoC) and  nucleotide counts for each genomic position of target genes for as
many samples as many BAM files are listed in the ".list" file.


The processing time depends on the dimension size of bam files and the number of
genes to investigate. If many genes need to be analyzed, we recommend using  
`buildInput` function before to launch  the App and to load the output file 
in Coverage Analysis page 

### Documentation 

To prepare input file and investigate coverage annotation the following 
**input parameters** must be specified in the sidebar of the **preprocessing** 
section:


- `Genome`: choose reference genome (hg19 or hg38) 

- ` Chromosome Notation ` : chromosome notation BAM. IN bam file `number` refers to chromosome 
                notation as 1, 2, ..., X,.M, `chr` refers to chr1, chr2, 
                ... chrX, chrM.

-  ` Load input file ` : loading txt file containing HGNC gene name(s)
                             (one per row) or a target bed 
                             
- `Choose the type of your input file` : selecting the type of input file                           

- ` Load bam file(s) list ` : loading path to file with ".list" extension containing 
                              absolute paths to BAM files (one per row)
                              
- `minimum mapping quality (MAPQ) `: default  value 1

- `minimum base quality (QUAL)` :  default value 1

Users can download `Statistical_Summary` report to obtain a coverage metrics per genes 
(`List of genes name`) or per amplicons (`Target Bed`) according to uploaded input 
file. 
The report summarizes following information: mean, median,
number of positions under 20x and percentage of position above 20x



