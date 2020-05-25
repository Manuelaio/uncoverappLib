Sys.setenv("R_TESTS" = "")

library(uncoverappLib)
library(testthat)
library(processx)

testthat::test_that(
  ##test if app is launching
  "app launches",{
    skip_on_cran()
    skip_on_travis()
   x <- processx::process$new(
      "R",
      c(
        "-e",
        "uncoverappLib::run.uncoverapp()"
     )
    )
    print(x)
    Sys.sleep(5)
    testthat::expect_true(print(x$is_alive()))
    x$kill()
  }
)

