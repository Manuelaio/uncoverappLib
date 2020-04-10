###make reactive database given reference genome
txdb= reactive({if (input$UCSC_Genome == "hg19"){
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene}
  else{
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene}
})

#################make reactive dataframe given input file

mydata <- reactiveVal()
observeEvent(input$file1, {

  tmp <- read.table(input$file1$datapath, header = input$header, stringsAsFactors = FALSE)
  colnames(tmp)[1:3]=c("chromosome","start","end")
  colnames(tmp)[-1:-3]=paste0("sample_", head(seq_along(tmp),-3))
  ## do whatever is needed to parse the data
  mydata(tmp)
})
example.file <- system.file(
  "extdata",
  "POLG.example.bed.gz",
  package = "uncoverappLib"
)
observeEvent(input$example_data, {
  polg= read.table(example.file)
  colnames(polg)= c("chromosome","start","end","sample_1")
  mydata(polg)
})
###Output of loaded file

#output$text<- #DT::renderDataTable({
 # renderDataTable({
#  validate(need(input$file1 != "", "Please, upload your file"))
 # mydata()})

###make reactive dataset given input choosed by users

mysample<-reactive({
  i= input$Sample
  #print(mydata())
  mydata() %>%
    dplyr:: select(chromosome, start, end,i) %>%
    dplyr::rename(value=i)
})

filtered_low<- reactive ({
  #print(mysample())
 # validate(
    #need(input$file1 != "", "Unrecognized data set: Please upload your file")
  #)
  mysample() %>%
    dplyr:: filter(chromosome == input$Chromosome,
                   value <= as.numeric(input$coverage_co))
})
filtered_high<- reactive ({
  mysample() %>%
    dplyr:: filter(chromosome == input$Chromosome,
                   value > as.numeric(input$coverage_co))

})

#######output dateset filtered

#output$text_cv <- DT::renderDataTable({
 # validate(
#    need(input$Gene_name != "" & input$Sample !="", "Please select all required input: Gene, Chromosome, Coverage threshold and Sample")
#  )
 # filtered_low()})

