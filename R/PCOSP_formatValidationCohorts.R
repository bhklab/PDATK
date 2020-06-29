#' Format the validation cohorts
#'
##TODO:: HEEWON - documentation
#'
#' @examples
#' data(PDACexpressionData)
#' cohortList <- formatValidationCohorts(PDACexpressionData)
#'
#' @param validationCohorts A named \code{list} of \code{data.frame}s containing
#'   independent patient cohorts?
#'
#' @param environment The \code{environment} into which variables will be
#'   assigned. This parameter is missing by default, which will cause the
#'   function to return a list instead of assigning variables to an environment.
#'
#' @return A named \code{list} containing the data matrix and group labels for
#'   each cohort in validationCohorts or assigns variables to the specified
#'   environment if one \code{envronment} parameter is included.
#'
#' @details We recommend using the \code{environment()} function to set the
#'   environment parameter. This will assign the matrix and group objects to
#'   your current environment as names of the cohorts with '_mat' and '_grp'
#'   appended to it.
#'
#' @export
formatValidationCohorts <- function(validationCohorts, environment) {

  if(!missing(environment)) {
    paste0("Assigning variables to the specified environment
           for cohorts: ", names(validationCohorts))
  } else {
    cohortList <- vector("list", length(validationCohorts))
  }

  for (i in seq_along(validationCohorts)) {
    cohort <- validationCohorts[[i]]
    cohort_mat <- convertCohortToMatrix(cohort)
    cohort_group_labels <- whichNotCensoredYearOne(cohort)
    cohort_group <- ifelse(as.numeric.factor(cohort$OS) >= 365,
                           1,
                           0)[cohort_group_labels]

    if (!missing(environment)) {
      assign(paste0(names(validationCohorts)[i], '_mat'),
             cohort_mat,
             envir=environment)
      assign(paste0(names(validationCohorts)[i], '_grp'),
             cohort_group,
             envir=environment)
      assign(paste0(names(validationCohorts)[i], '_grpIndex'),
             cohort_group,
             envir=environment)
    } else {
      cohortList[[i]] <- list("mat"=cohort_mat,
                              "grp"=cohort_group,
                              "grpIndex"=cohort_group_labels)
    }
  }

  if (missing(environment)) {
    names(cohortList) <- names(validationCohorts)
    return(cohortList)
  }
}
