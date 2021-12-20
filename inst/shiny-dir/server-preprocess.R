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
   if(input$type_file=="target_bed"){
      tmp_gene1<-read.table(input$gene1$datapath, stringsAsFactors = F)
      tmp_gene<- tmp_gene1[1:4]
      colnames(tmp_gene)<- c('chr', 'start', 'end', 'SYMBOL')
      #print(head(tmp_gene))
      ff<- tmp_gene }else{
         tmp_gene<- scan(file_gene$datapath, character(), quote = "")
         ff<- tmp_gene}
   #print(tmp_gene)
})

list_bam<- reactive({
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
   if(input$type_file=="target_bed"){
      errore<- data.frame(matrix(ncol = 2, nrow = 0))
   }else{
      gene_list1<- gene_list()
      my_gene_name<- OrganismDbi::select(org.Hs.eg.db, key= gene_list1, columns="ENTREZID",
                                         keytype="ALIAS")
      #print(my_gene_name)
      ID<- my_gene_name$ENTREZID
      b<- c()
      for (i in ID){
         if( ! is.element(i, keys(all_gene(), "GENEID")))
            b[i]<- i
      }
      errore<- subset(my_gene_name, my_gene_name$ENTREZID %in% b)
      colnames(errore)[1]= " the following HGNC official gene names are
    unrecognizable, please choose a correct name and reload a file"
      return(errore)
      #print(head(errore,n=12))
   }
})


for_bed<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_bam()))
      return(NULL)
   if(nrow(no_entrID())!=0)
      return(no_entrID())
   if(input$type_file=="target_bed"){
      for_bed<- gene_list() }else{
         gene_list1<- gene_list()
         my_gene_name<- OrganismDbi::select(org.Hs.eg.db, key= gene_list1,
                                            columns="ENTREZID", keytype="ALIAS")
         ID<- my_gene_name$ENTREZID
         pre<- data.frame()
         for (i in ID){

            txid <- OrganismDbi::select(all_gene(), i, "TXNAME", "GENEID")[["TXNAME"]]
            cds <-OrganismDbi::exonsBy(all_gene(), by="tx", use.names=TRUE)
            coor<- as.data.frame(cds[names(cds) %in% txid])
            coordinate<- cbind(coor,i)
            colnames(coordinate)[11]<- "ENTREZID"
            cood<- dplyr::inner_join(coordinate, my_gene_name, by="ENTREZID")
            pre<- rbind(cood,pre)
            pre$start<- pre$start -10
            pre$end<- pre$end +10
           # print(pre)
         }
         #print(head(pre))
         for_bed<- data.frame()
         for (row in 1:nrow(pre)){
            sel_cc<- data.frame(as.character(pre$seqnames[row]),
                                pre$start[row], pre$end[row],pre$ALIAS[row],
                                stringsAsFactors = FALSE)
            for_bed<- rbind(for_bed, sel_cc)

         }
         colnames(for_bed)<- c('chr', 'start', 'end', 'SYMBOL')
         for_bed<- unique(for_bed)
         if(input$notation == "number"){
            for_bed$chr<- gsub("^.{0,3}", "", for_bed$chr, perl =TRUE)}
         return(for_bed)
      }
})


pileup_input<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_bam()))
      return(NULL)
   if(nrow(no_entrID())!=0)
      return(no_entrID())

   for_grange<-GenomicRanges::makeGRangesFromDataFrame(for_bed(),
                                                       keep.extra.columns = TRUE)
   #print(head(for_grange))
   param <- Rsamtools::ScanBamParam(which= for_grange)

   p_param <- Rsamtools::PileupParam(distinguish_nucleotides=TRUE,
                                     distinguish_strands=FALSE,
                                     min_mapq=as.numeric(input$MAPQ),
                                     min_base_quality=as.numeric(input$base_qual),
                                     min_nucleotide_depth=1,
                                     max_depth=150000)
   df= list()
   for (i in list_bam()){
      pileup_df= Rsamtools::pileup(i, scanBamParam=param, pileupParam=p_param)

      df=rlist::list.append(df, pileup_df)
   }
   lst1 <- lapply(df, function(x) transform(x[,-5]))
   lst2<- lapply(lst1, function(x) transform(x[!duplicated(x),]))

   riarrange.df <- function(list_df){
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

   pp<- Reduce(function(...) merge(...,by= c("seqnames", "pos", "end"), all=TRUE), lst3)
   pp[is.na(pp)] <- 0
   pp<- as.data.frame(pp)
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


name_sample<- reactive({
   c<-gsub(".*/","",list_bam())
   return(c)
})

stat_summ<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_bam()))
      return(NULL)
   if(nrow(no_entrID())!=0)
      return(no_entrID())
   ppinp<- as.data.frame(pileup_input())
   colnames(ppinp)[1:3]<-c("seqnames","start","end")
   n<- length(colnames(ppinp)[-1:-3])
   #c<-gsub(".*/","",list_bam())
   m<-rep(name_sample(), each=2)
   colnames(ppinp)[-1:-3]<-paste0(c("sample_","nucleotide_"), m[1:n])
   for_grange<-GenomicRanges::makeGRangesFromDataFrame(for_bed(),
                                                       keep.extra.columns = TRUE)

   for_range_pp=GenomicRanges::makeGRangesFromDataFrame(ppinp,
                                                        keep.extra.columns = TRUE)
   tp<- GenomicRanges::findOverlaps(query= for_range_pp,
                                    subject = for_grange, type="any",
                                    select = "all")
   sts_df <- data.frame(for_range_pp[queryHits(tp),], for_grange[subjectHits(tp),])

   statistiche<- sts_df[!duplicated(sts_df$start),]
   statistiche<- subset(statistiche,
                        select = -c(width, strand, seqnames.1,
                                    start.1, end.1, width.1, strand.1))

   colnames(statistiche)[1:3]<- c("chromosome","start","end")
   merge_g<- dplyr::full_join(for_bed(),statistiche, by="SYMBOL", all=T)
   #col.sub<- colnames(merge_g[, grepl("sample_" , names(merge_g))])
   col_name= colnames(merge_g)

   col.sub= col_name[grepl("sample_", col_name)]
  merge_g[col.sub] <- sapply(merge_g[col.sub],as.numeric)

   x<- list()
   for (i in col.sub){
      x_df<- merge_g %>%
      dplyr::select("chromosome","start.x","end.x", "SYMBOL",i)
      #print(x_df)
      colnames(x_df)[5]<- "value"
      x_df$sample<- paste0(i)
      x<- list.append(x, x_df)
   }
   statistical_operation= function(df){
      df %>%
         dplyr::group_by(SYMBOL,sample) %>%
         dplyr::summarize(Total= sum(value, na.rm=TRUE),
                          Mean = mean(value, na.rm=TRUE),
                          Median= median(value, na.rm=TRUE),
                          number_of_position_under_20x = sum(value < 20),
                          percentage_under_20x= (sum(value < 20)/sum(value, na.rm=TRUE))*100) %>%
         as.data.frame(.) %>%
         dplyr::mutate_if(is.numeric, round, 3)
   }

   out_r<- do.call(rbind, lapply(x, function(x) statistical_operation(x)))
   return(out_r)


})
