#' Class Union for PCOSP, RLSModel and RGAModel Types
#'
#' @description
#' Class union used for method dispatch without inheritance
#'
setClassUnion('PCOSP_or_RLS_or_RGA', c('PCOSP', 'RLSModel', 'RGAModel'))

#' Class Union for PCOSP and ClinicalModel Types
#'
#' @description
#' Class union used for method dispatch without inheritance
#'
setClassUnion("PCOSP_or_ClinicalModel", c('PCOSP', 'ClinicalModel'))
