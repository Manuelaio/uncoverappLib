
test_that("build_input", {

  gene.list<- system.file("extdata", "mygene.txt", package = "uncoverappLib")
  bam_example <- system.file("extdata", "example_POLG.bam", package = "uncoverappLib")
  cat(bam_example, file = "bam.list", sep = "\n")
  temp_dir=tempdir()
  buildInput(geneList= gene.list, genome= "hg19",
                    type_bam= "chr",bamList= "bam.list", outDir= temp_dir )
  files <- list.files(temp_dir, full.names=TRUE)
  expect_true(file.exists(files))


  })
