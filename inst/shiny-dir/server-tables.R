#Tables
coord= eventReactive(input$button1,{
  disable("button1")
  shinyjs::show("text1")
  Sys.sleep(0.1)
  my_gene_name=OrganismDbi::select(org.Hs.eg.db,
                                   key= input$Gene_name,
                                   columns=c("ENTREZID","GENENAME", "ENSEMBL"),
                                   keytype="ALIAS")
  ID=my_gene_name$ENTREZID
  if (is.null(ID))
    return(NULL)
  all_gene= data.frame(genes(txdb()))
  pre= do.call(rbind, lapply(ID, function(x) data.frame(
    subset(all_gene, all_gene$gene_id == x), stringsAsFactors = FALSE)))
  colnames(pre)[6]= 'ENTREZID'
  info= merge(pre, my_gene_name)
  shinyjs::enable("button1")
  shinyjs::hide("text1")
  return(info)
})


observeEvent(input$button1, {
  output$ccg <-DT::renderDataTable({
    progress <- shiny::Progress$new()
    progress$set(message = "Running", detail = 'This may take a while')
    on.exit(progress$close())
    validate(
      need(input$Gene_name != "",
        "Unrecognized gene name: Please select HGNC gene name \n Click apply"))
    HGNC_org <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    validate(
      need(HGNC_org[HGNC_org %in% input$Gene_name],
           "incorrect HGNC Gene symbol"))
    coord()
    #print(coord())
  })
})

observeEvent(input$remove,{
  shinyjs::hide(coord())
  output$ccg<- NULL
})


exon_gp<-eventReactive(input$button1,{
  require(rtracklayer)
  ucsc <- browserSession()
  gname =input$Gene_name
  if (is.null(gname))
    return(NULL)
  else if (gname== "")
    return(NULL)
  if (input$UCSC_Genome == "hg19"){
    edb = EnsDb.Hsapiens.v75}
  else{
    edb= EnsDb.Hsapiens.v86}
  eid <- OrganismDbi::select(org.Hs.eg.db,gname,
                             "ENTREZID", "SYMBOL")[["ENTREZID"]]
  txid <- OrganismDbi::select(txdb(), eid, "TXNAME", "GENEID")[["TXNAME"]]
  cds <- cdsBy(txdb(), by="tx", use.names=TRUE)
  exoncds <- cds[names(cds) %in% txid]
  exon_Id<- as.data.frame(exoncds)
  exon_table=exon_Id  %>%
    dplyr::select(1:6,8,10)
  #exon_table=exon_Id[c(1:6,8,10)]
  colnames(exon_table)=c("number_of_transcript",
                         "type_of_transcript", "chrom", "start","end",
                         "length_of_exon", "cds_id", "exon_rank")
  return(exon_table)
  print(head(exon_table))
})

observeEvent(input$button1, {
  output$exon_pos<- DT::renderDataTable({
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress",
                 detail = 'This may take a while', value = 0)
    Sys.sleep(0.1)

    HGNC_org <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    validate(
      need(HGNC_org[HGNC_org %in% input$Gene_name],
           "incorrect HGNC Gene symbol"))
    exon_gp()
  })
})

observeEvent(input$remove,{
  output$exon_pos<- DT::renderDataTable({
    shinyjs::js$reset(exon_gp()) })
})


tryObserve <- function(x) {
  x <- substitute(x)
  env <- parent.frame()
  observe({
    tryCatch(eval(x, env),
             error = function(e) {
               #showNotification(paste("Error: ", e$message), type = "error")
             })
  })
}

tryObserve({
  if (is.null(input$exon_number) )
    return()
  message("load inputs")
  #observe({
  # if (is.null(input$exon_number))
  #  return(NULL)
  if (is.null(exon_gp()))
    return(NULL)
  exon_df= as.data.frame(exon_gp())
  one=subset(exon_df,
             exon_df$exon_rank == input$exon_number &
               exon_df$number_of_transcript == input$transcript_id)
  print(one)
  start_exon= as.numeric(one$start)
  end_exon= as.numeric(one$end)
  chr= as.numeric(gsub("\\D", "", one$chrom, perl =TRUE))
  trascript_name= as.character(one$type_of_transcript)
  updateTextInput(session, "id_t", value = trascript_name)
  updateTextInput(session, "Start_genomicPosition", value = start_exon)
  updateTextInput(session, "end_genomicPosition", value = end_exon)
  updateTextInput(session, "query_Database",
                  value = paste0(chr,":",start_exon,"-",end_exon))
})

#I define reactive START and END position

St<- reactive({
  sta= as.numeric(input$Start_genomicPosition)
  return(sta)
})

En<- reactive({
  en= as.numeric(input$end_genomicPosition)
  return(en)
})

id<- reactive({
  id_ucsc= as.character(input$id_t)
  return(id_ucsc)
})

coordinate<- reactive({
  coo=paste0(chr,":",start_exon,"-",end_exon)
  return(coo)
})
