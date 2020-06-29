#' PurIST
#'
#' Purity Independent Subtyping of Tumors
#'
#' @param dat \code{vector} A single sample gene expression profile 
#'     with gene symbols in names
#' 
#' @return probability
#'
#' @examples
#' probabibity <- PurIST(dat)
#'
#' @export
#' @importFrom gtools inv.logit
PurIST <- function(dat)
{
    biomarkers = c("GPR87", "REG4", 
                   "KRT6A", "ANXA10",
                   "BCAR3", "GATA6",
                   "PTGES", "CLDN18",
                   "ITGA3", "LGALS4",
                   "C16orf74", "DDC", 
                   "S100A2", "SLC40A1", 
                   "KRT5", "CLRN3")

    has_genes = is.element(biomarkers, names(dat))
    if (!all(has_genes)) {
        stop(sprintf("Please provide the following genes\n%s",
            paste(biomarkers, collapse = ", ")))
    }

    coefficient <- c(1.994, 2.031, 1.618, 0.922, 1.059, 0.929, 2.505, 0.485)
    intercept <- -6.815

    if (dat["GPR87"] > dat["REG4"])     { penal_1 <- 1 } else { penal_1 <- 0 }
    if (dat["KRT6A"] > dat["ANXA10"])   { penal_2 <- 1 } else { penal_2 <- 0 }
    if (dat["BCAR3"] > dat["GATA6"])    { penal_3 <- 1 } else { penal_3 <- 0 }
    if (dat["PTGES"] > dat["CLDN18"])   { penal_4 <- 1 } else { penal_4 <- 0 }
    if (dat["ITGA3"] > dat["LGALS4"])   { penal_5 <- 1 } else { penal_5 <- 0 }
    if (dat["C16orf74"] > dat["DDC"])   { penal_6 <- 1 } else { penal_6 <- 0 }
    if (dat["S100A2"] > dat["SLC40A1"]) { penal_7 <- 1 } else { penal_7 <- 0 }
    if (dat["KRT5"] > dat["CLRN3"])     { penal_8 <- 1 } else { penal_8 <- 0 }

    penal <- c(penal_1, penal_2, penal_3, penal_4, penal_5, penal_6, penal_7, penal_8)

    return( gtools::inv.logit(intercept + sum(coefficient * penal)) )
}
