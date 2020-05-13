# -------------------------------------------------------------------------
# Unit tests for PanCuRx utilty functions ---------------------------------
# -------------------------------------------------------------------------

# calcualteMedAbsDev ------------------------------------------------------
#data("cohortsCommonGenes")

cohorts <- readRDS('../data/cohortsCommonGenes.rds')       
test_that("rowMADdf returns the correct values", {
  expect_equal_to_reference(calcRowMADdf(cohorts$pcsi), file="utilities-calcRowMADdf.rds")
})
