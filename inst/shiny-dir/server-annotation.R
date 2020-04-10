##intersection between dbNSFP and filtered user-defined bed file
#ann.file<- system.file(
 # "extdata",
  #"sorted.bed.gz",
  #package = "uncoverappLib"
#)
#tbi<- system.file(
 # "extdata",
  #"sorted.bed.gz.tbi",
  #package = "uncoverappLib"
#)

intBED<- reactive({
  if (is.null(filtered_low()))
   return(NULL)
  bedA<- filtered_low()
  #file.name = (ann.file)
  #ANNOTATION FILEIN THE FOLDER OF SHINY SCRITP !!!
  #second and tirth columns are hg19 positions
  query.regions <- c(input$query_Database)
  if (is.null(query.regions))
    return(NULL)
  print(query.regions)
  result <- try({
    bedB <- tabix(query.regions, file.name, check.chr = FALSE)
  }, silent = TRUE)
  if ("try-error" %in% class(result)) {
    err_msg <- 'no coordinates recognized'}
  #print(err_msg)
  #print(head(bedB))
  colnames(bedB)<- c ('Chromo', 'start_hg19','end_hg19','REF','ALT',
                      'dbsnp','GENENAME', 'PROTEIN_ensembl', 'start_hg38',
                      'end_hg38','MutationAssessor','SIFT','Polyphen2',
                      'M_CAP','CADD_PHED','AF_gnomAD','ClinVar',
                      'clinvar_MedGen_id','clinvar_OMIM_id','HGVSc_VEP',
                      'HGVSp_VEP')
  #print(head(bedB))
  str(bedB)
  if (input$UCSC_Genome == "hg19"){
    bedB<- bedB %>%
      dplyr::rename(
        start="start_hg19",
        end="end_hg19"
      )
    #print(head(bedB))
  }
  else{
    bedB<- bedB[c("Chromo", "start_hg38","end_hg38","REF",
                  "ALT","dbsnp", "GENENAME", "PROTEIN_ensembl",
                  "start_hg19","end_hg19",'MutationAssessor','SIFT',
                  'Polyphen2','M_CAP','CADD_PHED','AF_gnomAD',
                  'ClinVar','clinvar_MedGen_id','clinvar_OMIM_id',
                  'HGVSc_VEP', 'HGVSp_VEP')]
    bedB<-bedB %>%
      dplyr::rename(
        start="start_hg38",
        end="end_hg38")
    #print(head(bedB))
    }

  for (i in bedB[1]){
    Chromosome<-paste("chr", i, sep="")
    bedB<- cbind(Chromosome, bedB)
    bedB[,2]<-NULL
  }
  #print(head(bedB))
  bedB$Chromosome= as.character(bedB$Chromosome)
  bedB$AF_gnomAD= suppressWarnings(as.numeric(bedB$AF_gnomAD))
  bedB$CADD_PHED= suppressWarnings(as.numeric(bedB$CADD_PHED))
  intersectBedFiles.GR <- function(bed1,bed2) {
    require(GenomicRanges)
    require(IRanges)
    bed1 <- makeGRangesFromDataFrame(bedA,ignore.strand = TRUE,
                                     keep.extra.columns = TRUE)
    bed2 <- makeGRangesFromDataFrame(bedB,ignore.strand = TRUE,
                                     keep.extra.columns = TRUE)
    tp<- findOverlaps(query = bed1, subject = bed2, type="any")
    intersect_df = data.frame(bed1[queryHits(tp),], bed2[subjectHits(tp),])
    return(intersect_df)
  }
  intersect_df<- intersectBedFiles.GR(bedA, bedB)
  return(intersect_df)
  #print(head(intersect_df))
})

#make reactive dataframe

condform_table<- reactive ({
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "table construction in progress",
               detail = 'This may take a while', value = 0)
  Sys.sleep(0.1)
  library(condformat)
  if (is.null(intBED()))
    return(NULL)
  condformat(intBED()) %>%
    rule_fill_discrete(ClinVar, expression= ClinVar !=".",
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(CADD_PHED, expression= CADD_PHED > 20,
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(MutationAssessor, expression=  MutationAssessor =='H',
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(M_CAP, expression=  M_CAP =='D',
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(AF_gnomAD, expression=
                         ifelse(is.na(AF_gnomAD) | AF_gnomAD < 0.5,
                                'TRUE','FALSE') ,
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(c(start, end),
                       expression = grepl("H|M", MutationAssessor) &
                         ClinVar !="." &  AF_gnomAD < 0.5 ,
                       colours= c("TRUE"= "yellow", "FALSE"= ""))%>%
    rule_css(c(start, end),
             expression = ifelse(grepl("H|M", MutationAssessor) &
                                   ClinVar !="." &  AF_gnomAD < 0.5,
                                 "red", "green"),
             css_field = "color")
})

