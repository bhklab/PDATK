setClassUnion('MAE_or_SE', c('MultiAssayExperiment', 'SummarizedExperiment'))

#' @title An `S4Model` Containing Molecular Data to be Consensus Clustered
#' 
#' @description
#' An `S4` wrapper class providing an interface to the `ConsensusClusterPlus`
#'   function from the package with the same name.
#' 
#' @inherit .S4Model
#' 
#' @md
#' @export
.ConsensusMetaclusteringModel <- setClass('ConsensusMetaclusteringModel',
    slots=c(trainData='MAE_or_SE'),
    contains='S4Model')

#' Constructor for a `ConsensusClusterModel` Object.
#' 
#' @param trainData A `MultiAssayExperiment` or `SummarizedExperiment` containing
#'   molecular data to be used for consensus clustering with 
#'   `ConsensusClusterPlus::ConsensusClusterPlus`
#' @param ... Force subsequent parameters to be named. Not used.
#' @param randomSeed An `integer` randomSeed that will be passed to the
#'   randomSeed parameter of the `ConsensusClusterPlus::ConsensusClusterPlus`
#'   function when training the model.
#' 
#' @return A `ConsensusMetaclusteringModel` object containing the training
#'   data and ready to be trained.
#' 
#' @seealso [`ConsensusClusterPlus::ConsensusClusterPlus`]
#' 
#' @aliases ConClustModel
#' 
#' @md
#' @export
ConsensusMetaclusteringModel <- function(trainData, ..., randomSeed) {
    funContext <- .context(1)

    if (missing(randomSeed)) stop(.errorMsg(funContext, 'No random seed was ',
        'specied for your model. Please include the value used for set.seed ',
        'when training this model! This ensures other can reproduce your ',
        'results.'))

    if (!is(trainData, 'MAE_or_SE')) stop(.errorMsg(funContext, 'The trainData',
        'argument is a ', class(trainData), ' object. It must be either a ',
        'MultiAssayExperiment or SummarizedExperiment'))

    return(.ConsensusMetaclusteringModel(
        trainData=trainData, 
        modelParams=SimpleList(randomSeed=randomSeed),
        models=SimpleList(),
        validationStats=data.table(),
        validationData=SimpleList(),
        metadata=list()
        ))
}
#' @export
ConMetaclustModel <- ConsensusMetaclusteringModel