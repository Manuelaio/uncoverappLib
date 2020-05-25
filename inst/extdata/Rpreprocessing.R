#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
output= args[3]
#setwd(output)
#dir.create('teMpFoldeR', showWarnings= FALSE)

mylist=scan(args[1], character(), quote = "")
#mylist=scan("./mylist.txt", character(), quote = "")
require(OrganismDbi)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Hs.eg.db)
require(Rsamtools)
#require(dplyr)
my_gene_name=OrganismDbi::select(org.Hs.eg.db, key= mylist, columns="ENTREZID",
                    keytype="ALIAS")


for (i in my_gene_name$ENTREZID){
    if (is.na(i)){
        s=print(subset(my_gene_name$ALIAS,is.na(my_gene_name$ENTREZID)))
        write.table(s, file='./preprocessing_log.txt', quote= FALSE,
                    row.names = FALSE, col.names = FALSE)
        stop('ERROR: unrecognized GENE NAME. Please choose a
         correct HGNC official gene names.')}
    else {return= NULL}
}

ID=my_gene_name$ENTREZID
genome= args[2]

if (genome == "hg19"){
    all_gene= TxDb.Hsapiens.UCSC.hg19.knownGene}else{
        all_gene= TxDb.Hsapiens.UCSC.hg38.knownGene}

b=c()
for (i in ID){
    if( ! is.element(i, keys(all_gene, "GENEID")))
        b[i]= i
}
no_entrID= subset(my_gene_name, my_gene_name$ENTREZID %in% b)

if (nrow(no_entrID)!=0){
    write.table(no_entrID, file='./preprocessing_log1.txt',
                quote= FALSE, row.names = FALSE, col.names = FALSE)
    stop('ERROR: unrecognized GENE NAME. Please remove genes
       stored in preprocessing_log1.txt')}


pre= data.frame()
for (i in ID){
    txid <- select(all_gene, i, "TXNAME", "GENEID")[["TXNAME"]]
    cds <-exonsBy(all_gene, by="tx", use.names=TRUE)
    coordinate= as.data.frame(cds[names(cds) %in% txid])
    pre= rbind(coordinate,pre)
    pre$start= pre$start -10
    pre$end= pre$end +10
}


for_bed=data.frame()
for (row in 1:nrow(pre)){
    sel_cc= data.frame(as.character(pre$seqnames[row]),
                       pre$start[row], pre$end[row], stringsAsFactors = FALSE)
    for_bed= rbind(for_bed, sel_cc)
}
colnames(for_bed)=c('chr', 'start', 'end')
for_bed= unique(for_bed)

type_bam= args[4]

if(type_bam == "number"){
    for_bed$chr= gsub("^.{0,3}", "", for_bed$chr, perl =TRUE)}

for_grange=makeGRangesFromDataFrame(for_bed, keep.extra.columns = TRUE)

param <- ScanBamParam(which= for_grange)

p_param <- PileupParam(distinguish_nucleotides=TRUE,
                       distinguish_strands=FALSE,
                       min_mapq=1,
                       min_base_quality=1,
                       min_nucleotide_depth=0,
                       max_depth=15000)
list= scan(args[5], character(), quote= "")
library(rlist)
df= list()
for (i in list){
    pileup_df= pileup(i, scanBamParam=param, pileupParam=p_param)

    df=rlist::list.append(df, pileup_df)
}
lst1 <- lapply(df, function(x) transform(x[,-5]))
lst2<- lapply(lst1, function(x) transform(x[!duplicated(x),]))

riarrange.df = function(list_df){
  require(dplyr)
  list_df %>%
    mutate(end= pos) %>%
    group_by(seqnames,pos,end,nucleotide) %>%
    summarise(count=sum(count)) %>%
    arrange(nucleotide) %>%
    summarise(value= paste(sum(count)),
              counts=paste(nucleotide, count, sep=':',collapse=';'))
}

lst3<- lapply(lst2, riarrange.df)

pp=Reduce(function(...) merge(...,by= c("seqnames", "pos", "end"), all=TRUE), lst3)


pp[is.na(pp)] <- 0
pp=as.data.frame(pp)
print(str(pp))
if(type_bam == "number"){
for (i in pp[1]){
  Chromosome<-paste("chr", i, sep="")
  pp<- cbind(Chromosome, pp)
  pp[,2]<-NULL
}
}
write.table(x=pp, file =paste0(format(Sys.time(), "%a_%b_%d_%X_%Y"),'.bed'),
                      quote=FALSE, sep="\t", eol = "\n", row.names = FALSE,
                   col.names = FALSE)
