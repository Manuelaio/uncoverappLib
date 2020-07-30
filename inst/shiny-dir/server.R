server <- function (input, output, session){
  options(shiny.maxRequestSize=30*1024^2)


  script1 <- system.file(
    "extdata",
    "Rpreprocessing.R",
    package = "uncoverappLib")


  #attach static scripts


  output$dependence = downloadHandler(filename="Rpreprocessing.R",
                                      content=function(file){
                                        file.copy(script1,file)
                                      })
  #source script to load dataset or example file
  source('server-preprocess.R', local= TRUE)
  output$input1<- renderDataTable({
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress",
                 detail = 'This may take a while', value = 0)
    Sys.sleep(0.1)
    validate(
      need(try(!is.null(pileup_input())), "please upload a file with HGNC
      gene names and absolute path(s) to BAM file"))
    pileup_input()
  })

  source('server-reactiveDF.R', local= TRUE)

  output$text<- #DT::renderDataTable({
    DT::renderDataTable({
      #validate(need(input$file1 != "", "Please, upload your file"))
      validate(need(ncol(mydata()) != "0", "Please, upload your file"))
      mydata()})

  output$text_cv <- DT::renderDataTable({
    validate(need(input$Gene_name != "" & input$Sample !="",
                  "Please select all required input: Gene, Chromosome,
                  Coverage threshold and Sample"))
    filtered_low()})

  #source script to reactive gene and exon table
  source('server-tables.R', local= TRUE)

  #source script to plot all gene coverage obteined by gviz
  #plot and to table low coverage position in each exon
  source('server-plots.R', local=TRUE)

  output$all_gene<- renderPlot({
    validate(
      need(ncol(mydata()) != "0", "Unrecognized data set: Please
           upload your file"))  # check error message
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Please wait a few minutes: Making plot",
                 detail = 'This may take a while', value = 0)
    for (i in 1:40) {
      progress$set(message = "Please wait a few minutes: Making plot",
                   detail = 'This may take a while', value = i)
      Sys.sleep(1.0)
    }
    Sys.sleep(0.1)
    p1()
  })

  output$df.l<- DT::renderDataTable({
    table1()})

  #plot with sequence reference

  output$sequence<- renderPlot({
    validate(
      need(ncol(mydata()) != "0",
           "Unrecognized data set: Please load your file"))
    validate(
      need(input$Start_genomicPosition < input$end_genomicPosition,
           "Please selct the right genomic position: end position is
           lower than start "))
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    for (i in 1:40) {
      progress$set(message = "Please wait a few minutes: Making plot",
                   detail = 'This may take a while', value = i)
      Sys.sleep(1.0)}
    on.exit(progress$close())
    Sys.sleep(0.1)
    p3()
  })

  #source script to dbNSFP annotation
  source('server-annotation.R', local=TRUE)
  output$tabExon<- DT::renderDataTable({
    validate(
      need(input$Gene_name != "" & input$Sample != "",
           "Please select all required input: Gene, Chromosome,
           Coverage threshold and Sample")
    )
    validate(
      need(try(!is.null(intBED())),'Unrecognized coordinates:
      Please change exon number input and be sure that input box
      "Region coordinates" is filled.
           Please check your R environment has annotation file loaded.'))
    #print(try(is.null(intBED())))
    if (is.null(intBED()))
      return(NULL)
    a= intBED() %>%
      dplyr::select(MutationAssessor,  ClinVar)
    df=data.frame(length(which(a$ClinVar!='.')),
                  length(which(grepl("H|M", a$MutationAssessor))))
    colnames(df)= c('Low coverage positions in ClinVar',
                    'Low coverage positions with High or Medium impact')
    return(df)

  })

  #table output wiht Condformat pkg

  output$uncover_position<- condformat::renderCondformat({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress", value = 0)
    Sys.sleep(0.1)
    validate(
      need(input$Gene_name != "" & input$Sample != "",
           "Please select all required input: Gene, Chromosome,
           Coverage threshold and Sample"))
    validate(
      need(try(!is.null(intBED())),'Unrecognized coordinates:
      Please change exon number input and be sure that input box
      "Region coordinates" is filled.
           Please check your R environment has annotation file loaded.')
    )
    if(is.null(condform_table()))
      return(NULL)
    #print(head(condform_table()))
    condform_table() %>%
     dplyr::select(seqnames, start, end, coverage, counts, REF, ALT,
                   dbsnp, GENENAME, PROTEIN_ensembl,  MutationAssessor,SIFT,
                    Polyphen2,M_CAP,CADD_PHED,AF_gnomAD,
                    ClinVar,clinvar_MedGen_id,clinvar_OMIM_id,HGVSc_VEP,
                  HGVSp_VEP)
    })

  #download data wiht conditionalFormatting

  output$downloadData <- downloadHandler(
    filename = function() {
      paste('download',input$file1,Sys.Date(), '.xlsx', sep='')
    },
    content = function(file){
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Sheet1")
      negStyle <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
      posStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      highlighted<-openxlsx::createStyle (fgFill = "yellow")
      hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF",
                        fontSize=12,
                        fontName="Arial Narrow", fgFill = "#4F80BD")
      openxlsx::writeData(wb, "Sheet1",condform_table(), headerStyle = hs)
      openxlsx::conditionalFormatting(wb, "Sheet1",cols=19,
                            rows=(1:nrow(condform_table())+1),rule='$S2=="H"',
                            style =negStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1",cols=19,
                            rows=(1:nrow(condform_table())+1),
                            rule='$S2!="H"', style =posStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1",cols=22,
                            rows=(1:nrow(condform_table())+1),
                            rule='$V2=="D"', style =negStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1",cols=22,
                            rows=(1:nrow(condform_table())+1),
                            rule='$V2!="D"', style =posStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1", cols=23,
                            rows=(1:nrow(condform_table())+1),
                            rule=">=20", style = negStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1", cols=23,
                            rows=(1:nrow(condform_table())+1),
                            rule="<20", style = posStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1", cols=24,
                            rows=(1:nrow(condform_table())+1),
                            rule="<=0.05", style = negStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1", cols=25,
                            rows=(1:nrow(condform_table())+1),
                            rule='!="."', style = negStyle)
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )


  # source code for calculator maxAF

  source('server-maxAF.R', local=TRUE)
  output$maxAF <- renderText({signif(data()[[1]],3)})
  output$uncoverPosition <- condformat::renderCondformat ({
    validate(
      need(ncol(mydata()) != "0",
           "Unrecognized data set: Please load your file"),
      need(input$Gene_name != "" &
             input$Sample != "", "Please select all required input: Gene,
           Chromosome, Coverage threshold and Sample")
    )
    validate(
      need(try(!is.null(uncover_maxaf())),
           'Unrecognized coordinates: Please change exon number input and
           be sure that input box "Region coordinates" is filled')
    )
    uncover_maxaf()
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress", value = 0)
    Sys.sleep(0.1)

    uncover_maxaf() %>%
      dplyr::select(seqnames, start, end, coverage, counts,REF,
                    ALT, dbsnp, GENENAME, PROTEIN_ensembl, MutationAssessor,
                    SIFT,Polyphen2,M_CAP,CADD_PHED,AF_gnomAD,
                    ClinVar,clinvar_MedGen_id,clinvar_OMIM_id,
                    HGVSc_VEP,HGVSp_VEP)
  })

  ###download excel with maxAF###

  output$download_maxAF <- downloadHandler(
    filename = function() {
      paste('download_uncover_maxAF',Sys.Date(), '.xlsx', sep='')
    },
    content = function(file){
      wb1 <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb1, "Sheet2")
      negStyle <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
      posStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      highlighted<-openxlsx::createStyle (fgFill = "yellow")
      hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF",
                        fontSize=12,
                        fontName="Arial Narrow", fgFill = "#4F80BD")
      openxlsx::writeData(wb1, "Sheet2",uncover_maxaf(), headerStyle = hs)
      openxlsx::conditionalFormatting(wb1, "Sheet2",cols=19,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='$S2=="H"', style =negStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2",cols=19,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='$S2!="H"', style =posStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2",cols=22,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='$V2=="D"', style =negStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2",cols=22,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='$V2!="D"', style =posStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2", cols=23,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule=">=20", style = negStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2", cols=23,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule="<20", style = posStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2", cols=24,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule="<=0.05", style = negStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2", cols=25,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='!="."', style = negStyle)
      openxlsx::saveWorkbook(wb1, file, overwrite = TRUE)
      #condformat2excel(condform_table(), file, sheet_name = "Sheet1",
      #                overwrite_wb = F, overwrite_sheet = T)
    }
  )

  #source script for binomial distribution given genomic position
  source('server-binomial.R', local=TRUE)
  output$bd<- renderPlot({
    library(stats)
    library(stats)
    p=seq(0,1,by=0.0001)

    ic<-qbinom(p = c(0.025, 0.975), size = df_subset(), prob =input$p)
    barplot(dbinom(x=1:df_subset(), size=df_subset(), p= input$p) ,
            main = paste(ic))
  })

  output$pbinom<- renderPlot({
    library(stats)
    #c=input$num_all
    p=seq(0,1,by=0.0001)
    Fx=pbinom(q=1:df_subset(), size=df_subset(), prob = input$p)
    npro=Fx[as.numeric(input$num_all)]
    leg= (1-npro)*100

    #leggend= (1-Fx[input$num_all]) *100
    #print(leggend)
    if (is.na(npro)){
      txt= "under threshold"}
    else{
      txt= paste("the probability of detecting more of",
                 input$num_all, "reads wiht variant alleles would be:",
                 leg, "%")
    }
    plot(1:df_subset(), Fx, type='h',xlab= "number of trials",
         main= txt)
    abline(v=input$num_all, col= "red")
  })

  output$ci<- renderText({
    #return(df_subset())
    ci_print<-qbinom(p = c(0.025, 0.975), size = df_subset(), prob =input$p)
    ci_in=c(ci_print[1]:ci_print[2])
    number= ci_print[1]
    number2=ci_print[2]
    thr= as.numeric(df_subset())
    # if (thr < input$coverage_co) {
    if (number2 < input$num_all) {
      print(paste("<span style=\"color:red\">according to the binomial probability model,
                  there is 95% probability to observe from </span>", number,
                  "<span style=\"color:red\">to</span>",
                  number2,
                  "<span style=\"color:red\"> variant-supporting reads .
                  </span>"))

    }else{
      print(paste("<span style=\"color:blue\">according to the binomial probability model,
                  there is 95% probability to observe from </span>", number,
                  "<span style=\"color:blue\">to</span>",
                  number2,
                  "<span style=\"color:blue\"> variant-supporting reads
                  .</span>"))
    }
  })
  observe({
    if (input$close > 0) stopApp()
  })
}
