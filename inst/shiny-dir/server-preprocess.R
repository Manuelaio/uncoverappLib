require(OrganismDbi)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Hs.eg.db)
require(rlist)
require(Rsamtools)

ff=reactiveValues()
gene_list <- reactive ({
  file_gene<-input$gene1
  if (is.null(input$gene1)) {return(NULL)}
  tmp_gene<- scan(file_gene$datapath, character(), quote = "")
  ff= tmp_gene
  #print(ff)
})

list_bam= reactive({
  if (is.null(input$bam_list)) {return(NULL)}
    tmp_bam <- scan(input$bam_list$datapath, character(), quote = "")
    #list_bam= (tmp_bam)
    #print(tmp_bam)
})


all_gene<- reactive({ if (input$Genome == "hg19"){
  all_gene<- TxDb.Hsapiens.UCSC.hg19.knownGene}else{
    all_gene<-TxDb.Hsapiens.UCSC.hg38.knownGene}
})

no_entrID<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_bam()))
      return(NULL)
   gene_list1= gene_list()
   my_gene_name= OrganismDbi::select(org.Hs.eg.db, key= gene_list1, columns="ENTREZID",
                                     keytype="ALIAS")
   #print(my_gene_name)
   ID=my_gene_name$ENTREZID
   b=c()
   for (i in ID){
      if( ! is.element(i, keys(all_gene(), "GENEID")))
         b[i]= i
   }
   errore= subset(my_gene_name, my_gene_name$ENTREZID %in% b)
   colnames(errore)[1]= " the following HGNC official gene names are
   unrecognizable, please choose a correct name and reload a file"
   return(errore)

})



pileup_input<- reactive({
    if (is.null(gene_list()))
        return(NULL)
    if (is.null(list_bam()))
        return(NULL)
   if(nrow(no_entrID())!=0)
      return(no_entrID())

    gene_list1= gene_list()
    my_gene_name=OrganismDbi::select(org.Hs.eg.db, key= gene_list1,
                                       columns="ENTREZID", keytype="ALIAS")
    ID=my_gene_name$ENTREZID
   pre= data.frame()
   for (i in ID){
      txid <- OrganismDbi::select(all_gene(), i, "TXNAME", "GENEID")[["TXNAME"]]
      cds <-exonsBy(all_gene(), by="tx", use.names=TRUE)
      coordinate= as.data.frame(cds[names(cds) %in% txid])
      pre= rbind(coordinate,pre)
      pre$start= pre$start -10
      pre$end= pre$end +10
      #print(pre)
   }
   print(head(pre))
   for_bed=data.frame()
   for (row in 1:nrow(pre)){
     sel_cc= data.frame(as.character(pre$seqnames[row]),
                        pre$start[row], pre$end[row], stringsAsFactors = FALSE)
     for_bed= rbind(for_bed, sel_cc)
   }
   colnames(for_bed)=c('chr', 'start', 'end')
   for_bed= unique(for_bed)


   if(input$notation == "number"){
     for_bed$chr= gsub("^.{0,3}", "", for_bed$chr, perl =TRUE)}

   for_grange<-GenomicRanges::makeGRangesFromDataFrame(for_bed,
                                                       keep.extra.columns = TRUE)
   #print(for_grange)
   param <- Rsamtools::ScanBamParam(which= for_grange)

   p_param <- Rsamtools::PileupParam(distinguish_nucleotides=TRUE,
                                     distinguish_strands=FALSE,
                                     min_mapq=1,
                                     min_base_quality=1,
                                     min_nucleotide_depth=1,
                                     max_depth=150000)


   df= list()
   for (i in list_bam()){
     pileup_df= Rsamtools::pileup(i, scanBamParam=param, pileupParam=p_param)

     df=rlist::list.append(df, pileup_df)
   }
   lst1 <- lapply(df, function(x) transform(x[,-5]))
   lst2<- lapply(lst1, function(x) transform(x[!duplicated(x),]))

   riarrange.df = function(list_df){
     require(dplyr)
     list_df %>%
        dplyr::mutate(end= pos) %>%
        dplyr::group_by(seqnames,pos,end,nucleotide) %>%
        dplyr::summarise(count=sum(count)) %>%
        dplyr::arrange(nucleotide) %>%
        dplyr::summarise(value= as.numeric(paste(sum(count))),
                 counts=paste(nucleotide, count, sep=':',collapse=';'))
   }

   lst3<- lapply(lst2, riarrange.df)

   pp=Reduce(function(...) merge(...,by= c("seqnames", "pos", "end"), all=TRUE), lst3)
   pp[is.na(pp)] <- 0
   pp=as.data.frame(pp)
   if(input$notation == "number"){
     for (i in pp[1]){
       Chromosome<-paste("chr", i, sep="")
       pp<- cbind(Chromosome, pp)
       pp[,2]<-NULL
     }
   }
   colnames(pp)<-NULL
   return(pp)
})


