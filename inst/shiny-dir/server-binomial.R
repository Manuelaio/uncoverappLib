#####output binomial distribution#####
df_subset <- reactive({
  validate(
    need(input$start_gp != "", "Please, select genomic intervals"))
  mysample()
  a <- subset(mysample(),start == as.numeric(input$start_gp)&
                end == as.numeric(input$end_gp))
  cov_ex<- as.numeric(a$coverage)
  print(cov_ex)
})




