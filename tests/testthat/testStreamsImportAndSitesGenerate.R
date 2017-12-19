library(adssn)
context("Basic import and site generation functions")

test_that(
  "importStreams() is able to import a SpatialStreamNetwork with no sites",
  {
    middlefork04 <- system.file("lsndata/MiddleFork04.ssn", package = "SSN")
    test.streams <- importStreams(middlefork04)
    expect_equal(class(test.streams)[1], "SpatialStreamNetwork")
    expect_equal(nrow(test.streams@obspoints@SSNPoints[[1]]@point.data), 0)
  }
)

test_that(
  "generateSites() is able to create observed sites on an empty SpatialStreamNetwork",
  {
    middlefork04 <- system.file("lsndata/MiddleFork04.ssn", package = "SSN")
    test.streams <- importStreams(middlefork04)
    obs.design <- prd.design <- rep(3, 2)
    with.sites <- generateSites(test.streams, "areaPI", binomialDesign(obs.design), binomialDesign(prd.design))
    expect_equal(nrow(with.sites@obspoints@SSNPoints[[1]]@point.data), 6)
    expect_equal(nrow(with.sites@predpoints@SSNPoints[[1]]@point.data), 6)
  }        
)
