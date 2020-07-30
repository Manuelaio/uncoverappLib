#####output binomial distribution#####
df_subset <- reactive({
  validate(
    need(input$start_gp != "" , "Please, select genomic position"))
  mysample()
  validate(need(
    input$start_gp %in% mysample()$start,
    "Please, select a valide genomic position"))

  a <- subset(mysample(),start == as.numeric(input$start_gp) &
                end == as.numeric(input$start_gp))
  cov_ex<- as.numeric(a$coverage)
  print(cov_ex)
})


