#' KRASdep
#'
#' Calculate L-score
#'
#' @param dat \code{vector} A single sample gene expression profile
#'        with gene symbols in names
#'
#' @return L-score
#'
#' @examples
#' lScore <- calcLScore(dat)
#'
#' @export
#'
calcLScore <- function(dat)
{
    # UP: KRAS-dependent genes, N=105
    upGenes <- c("ADAM8", "ADRB2", "ANGPTL4", "ARNTL2", "CALM2", 
                 "CALU", "CAPZA1", "CCL20", "CD274", "CDCP1", 
                 "CLCF1", "CSNK1D", "CXCL1", "CXCL2", "CXCL3", 
                 "CXCL5", "CXCL8", "DENND2C", "DUSP1", "DUSP4", 
                 "DUSP5", "DUSP6", "EFNB1", "EGR1", "EHD1", 
                 "ELK3", "EREG", "FERMT1", "FOS", "FOXQ1", 
                 "G0S2", "GDF15", "GLTP", "HBEGF", "IER3", 
                 "IL13RA2", "IL1A", "IL1B", "ITGA2", "ITPR3", 
                 "KCNK1", "KCNN4", "KLF5", "KLF6", "LAMA3", 
                 "LDLR", "LHFPL2", "LIF", "MALL", "MAP1LC3B", 
                 "MAP7D1", "MAST4", "MMP14", "MXD1", "MYDGF", 
                 "NAMPT", "NAV3", "NDRG1", "NFKBIZ", "NIPAL1", 
                 "NT5E", "OXSR1", "PHLDA1", "PHLDA2", "PI3", 
                 "PIK3CD", "PIM1", "PLAUR", "PNMA2", "PPP1R15A", 
                 "PRNP", "PTGS2", "PTHLH", "PTPRE", "PTX3", 
                 "PVR", "S100A6", "SDC1", "SDC4", "SEMA4B", 
                 "SERPINB1", "SERPINB2", "SERPINB5", "SESN2", "SFN", 
                 "SLC16A3", "SLC2A14", "SLC2A3", "SLC9A1", "SPRY4", 
                 "TFPI2", "TGFA", "TIMP1", "TMEM45B", "TNFRSF10A", 
                 "TNFRSF10B", "TNFRSF12A", "TNS4", "TOR1AIP1", "TSC22D1", 
                 "TUBA4A", "UAP1", "UPP1", "VEGFA", "ZFP36")

    # DOWN: KRAS-independent genes, N=41
    downGenes <- c("ABCC5", "ACAP2", "ARMC8", "ATPAF1", "AUTS2",
                   "CCSAP", "CELSR2", "CEP57L1", "COQ7", "DRD4", 
                   "ENAH", "ERP44", "GREB1L", "HNRNPU", "HTATSF1", 
                   "ID4", "ITSN1", "KDM4C", "MIB1", "MRPS14", 
                   "MSI1", "MSI2", "NUP133", "OGN", "PARP1", 
                   "PIAS1", "RASL10B", "RTN3", "SEC63", "SH3GL2", 
                   "SMAD9", "STARD7", "SUGP1", "TBC1D24", "TMEFF1", 
                   "TTC28", "ZNF292", "ZNF441", "ZNF493", "ZNF669", 
                   "ZNF672")

    hasUpGenes <- is.element(upGenes, names(dat))
    hasDownGenes <- is.element(downGenes, names(dat))

    if (!all(hasUpGenes)) {
        stop(sprintf("Please provide the following genes: %s",
                      paste(upGenes[!hasUpGenes], collapse=", ")))
    } else if (!all(hasDownGenes)) {
        stop(sprintf("Please provide the following genes: %s", 
                      paste(downGenes[!hasUpGenes], collapse=", ")))
    } else {
        upGenesIdx <- which(names(dat) %in% upGenes)
        downGenesIdx <- which(names(dat) %in% downGenes)
    }

    return(mean(dat[upGenesIdx]) - mean(dat[downGenesIdx])) 
}
#' Calculate S-score
#'
#' @param dat \code{vector} A single sample gene expression profile
#'        with gene symbols in names
#'
#' @return S-score
#'
#' @examples
#' sScore <- calcSScore(dat)
#'
#' @export
#'
calcSScore <- function(dat)
{
    # UP: KRAS-dependent genes, N=79
    upGenes <- c("ACTR3C", "ANKRD22", "ATP2C2", "B3GNT3", "BSPRY", 
                 "C11orf52", "C1orf106", "C1orf116", "C1orf172", "C1orf210", 
                 "C1orf74", "C6orf141", "CDA", "CDH1", "CDS1", 
                 "CEACAM6", "CGNL1", "CLDN4", "CLDN7", "DENND1C", 
                 "DNAJA4", "DPP4", "DSC2", "EFNB2", "EHF", 
                 "ELF3", "EPCAM", "EPN3", "EPS8L1", "ERBB3", 
                 "ESRP1", "ESRP2", "F11R", "FRK", "GALNT3", 
                 "GPR115", "GRHL2", "HS3ST1", "INPP4B", "IRF6", 
                 "ITGB6", "KRTCAP3", "LAD1", "LAMA3", "MAL2", 
                 "MAOA", "MAPK13", "MPZL2", "MPZL3", "MYO1D", 
                 "OVOL2", "PCDH1", "PGM2L1", "PLEKHA7", "PPL", 
                 "PPP1R14C", "PROM2", "RAB11FIP4", "RAB17", "RAB25", 
                 "RALGPS2", "S100A14", "SCEL", "SCIN", "SCNN1A", 
                 "SDR16C5", "SPINT1", "ST14", "SYK", "SYTL5", 
                 "TGFA", "TIAF1", "TMEM154", "TMEM30B", "TMEM45B", 
                 "TMPRSS4", "TSPAN1", "TTC9", "VSIG1")

    # DOWN: KRAS-independent genes, N=40
    downGenes <- c("ABP1", "ACTA2", "ADRA2C", "ANTXR1", "ANXA6", 
                   "APLN", "DYRK3", "EML1", "FAR1", "HNRNPA2B1", 
                   "HSPA12A", "HTRA1", "IKBIP", "KCTD15", "LIX1L", 
                   "LOC541471", "MAGEE1", "MPHOSPH9", "MPPE1", "MSRB3", 
                   "NUDT11", "OSTM1", "PARVB", "PAX6", "PLCB4", 
                   "PPARGC1B", "PPP4R2", "RECK", "RHOT1", "RYK", 
                   "SLC1A3", "SLC47A1", "SMARCD3", "SRGN", "SYDE1", 
                   "SYNGR1", "TMEM237", "TUB", "TXNRD1", "WDR35")

    hasUpGenes <- is.element(upGenes, names(dat))
    hasDownGenes <- is.element(downGenes, names(dat))

    if (!all(hasUpGenes)) {
        stop(sprintf("Please provide the following genes: %s",
                      paste(upGenes[!hasUpGenes], collapse=", ")))
    } else if (!all(hasDownGenes)) {
        stop(sprintf("Please provide the following genes: %s", 
                      paste(downGenes[!hasUpGenes], collapse=", ")))
    } else {
        upGenesIdx <- which(names(dat) %in% upGenes)
        downGenesIdx <- which(names(dat) %in% downGenes)
    }

    return(mean(dat[upGenesIdx]) - mean(dat[downGenesIdx])) 
}
#' Transform raw TPM to log10(TPM)
#'
#' @param mat \code{matrix} A expression profile in raw TPM values
#'
#' @return Relative TPM values
#'
#' @examples
#' data <- transExp(mat)
#'
#' @export
#'
transExp <- function(mat)
{
    if (all(is.na(mat)) || ncol(mat) < 2 || nrow(mat) < 2) {
        stop("Please provide raw TPM values in matrix.")
    }
    # log10-transformation
    log10TPM <- log10(mat)

    # average log10(TPM) per gene
    avgLog10TPM <- apply(log10TPM, 1, mean)

    # get relative values (log10(TPM) - mean(log10(TPM)))
    return(log10TPM - avgLog10TPM)
}
