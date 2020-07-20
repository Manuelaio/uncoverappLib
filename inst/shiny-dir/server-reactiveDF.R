###make reactive database given reference genome
txdb= reactive({if (input$UCSC_Genome == "hg19"){
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene}
  else{
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene}
})

#################make reactive dataframe given input file

mydata <- reactiveVal()
observeEvent(input$file1, {

  tmp <- read.table(input$file1$datapath,
                    header = input$header, stringsAsFactors = FALSE)
  colnames(tmp)[1:3]=c("chromosome","start","end")
  n=  length(colnames(tmp)[-1:-3])
  a=rep(head(seq_along(tmp),-3), each=2)
  colnames(tmp)[-1:-3]=paste0(c("sample_","nucleotide_"), a[1:n])
  ## do whatever is needed to parse the data
  mydata(tmp)
})

observeEvent(input$pileup, {

  tmp_pileup <- pileup_input()
  colnames(tmp_pileup)[1:3]=c("chromosome","start","end")
  n=  length(colnames(tmp_pileup)[-1:-3])
  a=rep(head(seq_along(tmp_pileup),-3), each=2)
  colnames(tmp_pileup)[-1:-3]=paste0(c("sample_","nucleotide_"), a[1:n])
  ## do whatever is needed to parse the data
  mydata(tmp_pileup)
})

#example.file <- system.file(
 # "extdata",
  #"POLG.example.bed",
  #package = "uncoverappLib"
#)

#observeEvent(input$example_data, {
 # polg= read.table(example.file)
  #colnames(polg)= c("chromosome","start","end","sample_1","nucleotide_1" )
  #mydata(polg)
#})
###Output of loaded file

###make reactive dataset given input choosed by users

mysample<-reactive({
  if (is.null(mydata()))
    return(NULL)
  #i= input$Sample
  num= input$Sample
  i=paste0("sample_",num)
  nucleodites= paste0("nucleotide_",num)
  #print(mydata())
  mydata() %>%
    dplyr:: select(chromosome, start, end,i,nucleodites) %>%
    dplyr::rename(coverage=i) %>%
    dplyr::rename(counts= nucleodites)
})


filtered_low<- reactive ({
  if (is.null(mysample()))
    return(NULL)
  mysample() %>%
    dplyr::select(-c(counts)) %>%
    dplyr::filter(chromosome == Chromosome(),
                  coverage <= as.numeric(input$coverage_co))

})

filtered_high<- reactive ({
  if (is.null(mysample()))
    return(NULL)
  mysample() %>%
    dplyr::select(-c(counts)) %>%
    dplyr::filter(chromosome == Chromosome(),
                  coverage > as.numeric(input$coverage_co))
})
