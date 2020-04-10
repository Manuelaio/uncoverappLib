#Max calculation

data <- reactive({
  myPrev = 1/input$prev
  if(input$inh=="monoallelic"){
    myMaxAF = (1/2) * myPrev * input$hetA * input$hetG * (1/input$pen)
  }
  if(input$inh=="biallelic"){
    myMaxAF = sqrt(myPrev) * input$hetA * sqrt(input$hetG) * (1/sqrt(input$pen))
  }
  myMaxAC = qpois(p=as.numeric(input$CI),
                  lambda=(input$popSize)*(myMaxAF))
  return(list(myMaxAF,myMaxAC))
})

#output$maxAF <- renderText({signif(data()[[1]],3)})

uncover_maxaf <- reactive ({

  condformat(intBED()) %>%
    rule_fill_discrete(ClinVar, expression= ClinVar !=".",
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(MutationAssessor, expression= MutationAssessor =='H',
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(c(M_CAP, AF_gnomAD),
                       expression= M_CAP > 0.025 & AF_gnomAD <0.05,
                       colours= c("TRUE"= "red", "FALSE"= "green")) %>%
    rule_fill_discrete(c(start, end),
                       expression= MutationAssessor =='H' &
                         ClinVar !="." & M_CAP != "TRUE" &
                         AF_gnomAD < (signif(data()[[1]],3)) ,
                       colours= c("TRUE"= "yellow", "FALSE"= ""))%>%
    rule_css(c(start, end),
             expression= ifelse(grepl("H|M", MutationAssessor) &
                                  ClinVar !="." & M_CAP != "TRUE" &
                                  AF_gnomAD < (signif(data()[[1]],3)),
                                "red", "green"), #& CADD_PHED > 20
             css_field = "color")

})


