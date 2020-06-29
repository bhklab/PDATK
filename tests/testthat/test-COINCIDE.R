# -------------------------------------------------------------------------
# COINCIDE implementation tests ------------------------------------------
# -------------------------------------------------------------------------

## FIXME:: Do this on data subset!

# calcGeneWiseWeightedMADs ------------------------------------------------
# FIXME:: Load in necessary data
test_that("calcGeneWiseWeightedMADs returns correct values", {
  expect_equal_to_reference(calcGeneWiseWeightedMADs(cohortMADrankings), 
                      file='COINCIDE-calcGeneWiseWeightedMADs.rds')
})


# getTopGenes -------------------------------------------------------------
## FIXME:: Load in necessary data
test_that("getTopGenes returns the correct values", {
  expect_equal_to_reference(getTopGenes(cohortMADrankings, 20), 
                      file="COINCIDE-getTopGenes.rds")
})


# getMetaGenes ------------------------------------------------------------
test_that("getMetaGenes returns the correct values", {
  expect_equal_to_reference(getMetaGenes(cohortMADrankings, geneWiseWeightedMADdf, 
                                   n=20, m=2066), 
                      file="COINCIDE-getMetaGenes.rds")
})


# .subsetCohorts ----------------------------------------------------------
test_that(".subsetCohorts returns the correct values", {
  expect_equal_to_reference(.subsetCohorts(cohorts, genes=metaGenes), 
                      file="COINCIDE-dot-subsetCohorts.rds")
})

# .scaleGenewise ----------------------------------------------------------
## FIXME:: Add data
test_that(".scaleGenewise returns the correct values", {
  cohortSubset <- .subsetCohorts(cohorts, genes=metaGenes)[[1]]
  expect_equal_to_reference(.scaleGenewise(cohortSubset, center=TRUE, scale=FALSE), 
                      file="COINCIDE-dot-scaledGenewise.rds")
})


# preprocessCohorts -------------------------------------------------------
test_that("preprocessCohorts returns the correct values", {
  expect_equal_to_reference(preprocessCohorts(cohorts[1], metaGenes), 
                      file="COINCIDE-preprocessCohorts.rds")
})


# consensusClusterCohort --------------------------------------------------
test_that("consensusClusterCohort returns the correct values", {
  procCohort <- preprocessCohorts(cohorts[1], metaGenes)[[1]]
  expect_equal_to_reference(consensusClusterCohort(procCohort,
                                                   maxK=3,
                                                   distance="pearson",
                                                   method="hc"), 
                            file="COINCIDE-consensusClusterCohorts.rds")
})