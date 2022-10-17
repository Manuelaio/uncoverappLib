
#strack<- reactive({if (input$UCSC_Genome == "hg19"){
#  Gviz::SequenceTrack(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosome = Chromosome())}
#  else{
#    Gviz::SequenceTrack(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosome = Chromosome())}
#})

dtrack1<- reactive ({grcoverage <- filtered_low()
dt1<-Gviz::DataTrack(range= grcoverage,type = "histogram", fill.histogram =  "red", col.histogram="NA",
               genome = input$UCSC_Genome,name = "Seq. Depth")
})

dtrackHigh<- reactive ({ grcoverage_high <- filtered_high()
dt2<- Gviz::DataTrack(range = grcoverage_high,type = "histogram", fill.histogram = "dodgerblue",
                col.histogram="NA", genome = input$UCSC_Genome,name = "Seq. Depth")
})

itrack <- reactive ({
  i= Gviz::IdeogramTrack(genome =input$UCSC_Genome , chromosome = Chromosome())
})

grtrack  <-reactive({
  ggT<-Gviz::GeneRegionTrack (txdb(), chromosome  = Chromosome() ,  start  =  St() ,  end  =  En(),
                        showId  =  TRUE , transcriptAnnotation="symbol",
                        name  =  " Gene Annotation ")
  return(ggT)
})
#Prepration for all gene coverage plot
p1<- reactive ({
  start_gene= coord()$start[coord()$seqnames ==Chromosome()]
  if (length(start_gene) > 1)
    start_gene= start_gene[1]
  print(start_gene)
  end_gene= coord()$end[coord()$seqnames ==Chromosome()]
  if (length(end_gene) > 1)
    end_gene= end_gene[1]
  print(end_gene)

  ot <- Gviz::OverlayTrack(trackList = list(dtrack1(), dtrackHigh()))
  gtrack <- Gviz::GenomeAxisTrack()
  #ylims <- extendrange(range(c(values(dtrack1()), values(dtrackHigh()))))
  ylims <- grDevices::extendrange(range(c(dtrack1()@data), dtrackHigh()@data))
  gr_ex_track <- Gviz::GeneRegionTrack(txdb(),
                                       chromosome = Chromosome(),
                                       start = start_gene, end = end_gene,
                                 showId = TRUE,
                                 name = "Gene Annotation")
  p1= Gviz::plotTracks(list(itrack(), gtrack, ot, gr_ex_track),
                 from = start_gene,
                 to=end_gene,reverseStrand = TRUE,
                 ylim=ylims,type = "histogram",baseline=input$coverage_co,
                 col.baseline = "black",lwd.baseline=0.3,
                 extend.left = 0.5, extend.right = 200)

})

#table of uncovered exons

table1<- reactive({
  if (is.null(filtered_low()))
    return(NULL)
  f.low=filtered_low()
  f.low[,'new']='NA'
  f.low$new<- ifelse(sapply(f.low$start, function(p)
    any(exon_gp()$start <= p & exon_gp()$end >= p)), "YES", "out")
  l.coverage= as.data.frame(f.low[f.low$new=='YES',])
  validate(
    need(nrow(l.coverage) >0, "ALL EXONS ARE COVERED
         UNDER YOUR CHOOSE THRESHOLD"))

  x= l.coverage$start
  getValue3 <- function(x, data) {
    tmp <- data %>%
      dplyr::filter(start <= x, x <= end) %>%
      dplyr::filter(number_of_transcript == input$transcript_id)
    return(tmp$exon_rank)
  }

  a_exon=sapply(x, getValue3, data=exon_gp())
  exon=unlist(lapply(a_exon, function (x) ifelse(length (x) > 0, x, NA)))
  exon.df= as.data.frame(exon)
  if (is.null(exon.df))
    return(NULL)
  df.l= cbind(l.coverage, exon.df)
  t1=table(df.l$exon)
  df.t1= as.data.frame(t1)
  colnames(df.t1)= c('exon','uncovered positions')
  return(df.t1)
})

#output$df.l<- renderDataTable({
 # table1()
#})

#Preparation of exon zooom plot

p2<- eventReactive(input$button5, {
  gname =input$Gene_name
  if (is.null(gname))
    return(NULL)
  disable("button5")
  shinyjs::show("text1")
  Sys.sleep(0.1)
  gtrack <- Gviz::GenomeAxisTrack()
  ot <- Gviz::OverlayTrack(trackList = list(dtrack1(), dtrackHigh()))
  ylims <- grDevices::extendrange(range(c(dtrack1()@data), dtrackHigh()@data))
  one_trascript= grtrack()[grtrack()@range$transcript== id()]
  print(head(one_trascript))
  grtrack_symbol <- Gviz::GeneRegionTrack(one_trascript@range,
                                          chromosome = Chromosome(),
                                    start = St(),
                                    end = En(),
                                    showId = TRUE, exonAnnotation="exon",
                                    name = "Gene Annotation & Symbol")
  grtrack_range <- grtrack_symbol@range
  range_mapping <- OrganismDbi::select(Homo.sapiens,
                                       keys = mcols(grtrack_range)$symbol,
                                       keytype = "TXNAME",
                                       columns = c("ENTREZID", "SYMBOL"))
  library(stringr)
  new_symbols <- with(range_mapping,str_c(SYMBOL, " (", TXNAME, ")", sep = ""))
  symbol(grtrack_symbol) <- new_symbols
  shinyjs::enable("button5")
  shinyjs::hide("text1")
  p2=Gviz::plotTracks(list(itrack(), gtrack, ot, grtrack_symbol),
                      ylim=ylims,xlim=NULL,exonAnnotation="exon1",
                from = St(), to = En(), reverseStrand = TRUE,
                background.panel = "#FFFEDB", background.title = "darkblue",
                baseline=input$coverage_co,
                col.baseline = "black",lwd.baseline=0.3,
                extend.left = 0.5, extend.right = 100,
                main= paste0('exon',input$exon_number))

})

observeEvent(input$button5,{
  output$ens<- renderPlot({
    validate(
      need(ncol(mydata()) != "0", "Unrecognized data set:
           Please load your file"))
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    for (i in 1:40) {
      progress$set(message = "Please wait a few minutes:
                   Making plot", detail = 'This may take a while', value = i)
      Sys.sleep(1.0)
    }
    on.exit(progress$close())
    Sys.sleep(0.1)
    p2()
  })
})

observeEvent(input$remove5,{
  shinyjs::hide(p2())
  output$ens<- NULL
})

observeEvent(input$btn_run,{
  Sys.sleep(5)
})

####output plot with reference sequqnce, few bases##########

p3<- reactive({
  ot <- Gviz::OverlayTrack(trackList = list(dtrack1(), dtrackHigh()))
  gtrack <- Gviz::GenomeAxisTrack()
  ylims <- grDevices::extendrange(range(c(dtrack1()@data), dtrackHigh()@data))
  #strack<- SequenceTrack(homo, chromosome = Chromosome())
  #p3=Gviz::plotTracks(list(itrack(), gtrack, ot,grtrack(), strack()),
  p3=Gviz::plotTracks(list(itrack(), gtrack, ot,grtrack()),
                from = as.numeric(input$Start_genomicPosition),
                to=as.numeric(input$end_genomicPosition),
                reverseStrand = TRUE, cex = 0.8, ylim= ylims,
                baseline=input$coverage_co,
                col.baseline = "black",lwd.baseline=0.5,
                extend.right = 100, extend.left = 0.5)
})



