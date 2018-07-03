#' Deconvolve cell type mixing proportions from gene expression data. 
#' @param Y Expression matrix.
#'
#' (Required) Two-dimensional numeric. Must implement \code{as.matrix}.
#'
#' Each row contains expression measurements for a particular sample. Each columm contains the measurements of the same gene over all individuals. Can either contain just the mixture samples to be deconvolved or both the mixture samples and the reference samples. See \code{pure_samples} and \code{references} for more details.
#' @param references Cell-type reference expression matrix.
#'
#' (Optional) Two-dimensional numeric. Must implement \code{as.matrix}. Must have same number of columns as \code{Y}. Columns must correspond to columns of \code{Y}.
#'
#' Each row contains expression measurements for a reference profile of a particular cell type. Columns contain measurements of reference profiles of a gene. Optionally may merge this matrix with \code{Y} and use \code{pure_samples} to indicate which rows of \code{Y} are pure samples. If \code{pure_samples} is not specified \code{references} must be specified. In this case each row of \code{references} is assumed to be a distinct cell-type. If both \code{pure_samples} and \code{references} are specified then \code{pure_samples} specifies to which cell-type each row of \code{references} corresponds. 
#' @param pure_samples The pure sample indicies.
#'
#' (Optional) List of one-dimensional integer. Must implement \code{as.list}.
#'
#' The i-th element of the top-level list is a vector of indicies (rows of \code{Y} or \code{references}) that are pure samples of type i. If \code{references} is not specified then this argument identifies which rows of \code{Y} correspond to pure reference samples of which cell-types. If \code{references} is specified then this makes same idenficiation but for the \code{references} matrix instead. 
#' @param data_type Type of expression measurements.
#'
#' (Optional) One-dimensional string. 
#'
#' An optional string indicating the type of the expression measurements. This is used to set gamma to a pre-determined value based upon the data type. Valid values are for probe-level microarray as ``microarray-probe'', gene-level microarray as ``microarray-gene'' or rna-seq as ``rna-seq''. Alternatively can set \code{gamma} directly. 
#' @param n_markers Number of marker genes.
#'
#' (Optional) One-dimensional numeric.
#'
#' How many markers genes to use for deconvolution. Can either be a single integer, vector of integers (one for each cell type), or single or vector of percentages (numeric in 0 to 1). If a single integer then all cell types use that number of markers. If a vector then the i-th element determines how many marker genes are used for the i-th cell type. If single percentage (in 0 to 1) then that percentage of markers are used for all types. If vector of percentages then that percentage used for each type, respectively. If not specified then top 10\% of genes are used.
#' @param gamma Expression adjustment term.
#'
#' (Optional) One-dimensional positive numeric.
#'
#' If provided as a single positive number then that value will be used for \code{gamma} and over-ride the value of gamma chosen by the \code{data_type} argument. If neither \code{gamma} nor \code{data_type} are specified then \code{gamma} will be set to one. 
#' @param markers Marker gene indices.
#'
#' (Optional) List of one-dimensional integer.
#'
#' Top-level list should be same length as \code{pure_samples}, i.e. one element for each cell type. Each element of the top-level list is a vector of indicies (columns of \code{Y}) that will be considered markers of that particular type. If not supplied then \code{dtangle} finds markers internally using \code{find_markers}. Alternatively, one can supply the output of \code{find_markers} to the markers argument. 
#' @param marker_method Method used to rank marker genes.
#'
#' (Optional) One-dimensional string.
#' 
#' The method used to rank genes as markers. If not supplied defaults to ``ratio''. Only used if markers are not provided to argument ``markers''. Options are
#' \itemize{
#' \item{'ratio'}{ selects and ranks markers by the ratio of the mean expression of each gene in each cell type to the mean of that gene in all other cell types.}
#' \item{'regression '}{ selects and ranks markers by estimated regression coefficients in a series of regressions with single covariate that is indicator of each type.}
#' \item{'diff'}{ selects and ranks markers based upon the difference, for each cell type, between the median expression of a gene by each cell type and the median expression of that gene by the second most highly expressed cell type.}
#' \item{'p.value'}{ selects and ranks markers based upon the p-value of a t-test between the median expression of a gene by each cell type and the median expression of that gene by the second most highly expressed cell type.}
#' }
#' @param summary_fn What summary statistic to use when aggregating expression measurements.
#'
#' (Optional) Function that takes a one-dimensional vector of numeric and returns a single numeric.
#'
#' Defaults to the mean. Other good options include the median. 
#' @return List.
#' \itemize{
#' \item{'estimates'}{ a matrix estimated mixing proportions. One row for each sample, one column for each cell type.}
#' \item{'markers'}{ list of vectors of marker used for each cell type. Each element of list is vector of columns of \code{Y} used as a marker for the i-th cell type.}
#' \item{'n_markers'}{ vector of number of markers used for each cell type.}
#' \item{'gamma'}{ value of the sensitivity parameter gamma used by dtangle.}
#' }
#' @examples
#' truth = shen_orr_ex$annotation$mixture
#' pure_samples <- lapply(1:3, function(i) {
#'    which(truth[, i] == 1)
#' })
#' Y <- shen_orr_ex$data$log
#' n_markers = 20
#'
#' dtangle(Y, pure_samples = pure_samples,
#' n_markers=n_markers,data_type='microarray-gene',marker_method = 'ratio')
#'
#' n_markers = c(10,11,12)
#' dtangle(Y, pure_samples=pure_samples,
#' n_markers=n_markers,gamma=.8,marker_method = 'regression')
#' @seealso \code{\link{find_markers}}
#' @export
dtangle <- function(Y, references = NULL, pure_samples = NULL, n_markers = NULL, 
    data_type = NULL, gamma = NULL, markers = NULL, marker_method = "ratio", summary_fn = mean) {
    
    stopifnot(all(n_markers > 0))
    stopifnot(!is.null(c(references, pure_samples)))
    
    cmbd <- combine_Y_refs(Y, references, pure_samples)
    Y <- cmbd$Y
    pure_samples <- cmbd$pure_samples
    
    prc <- process_markers(Y, pure_samples, n_markers, data_type, gamma, markers, 
        marker_method)
    n_markers <- prc$n_markers
    mrkrs <- prc$mrkrs
    gamma <- prc$gamma
    
    baseline <- baseline_exprs(Y, pure_samples, mrkrs, summary_fn = summary_fn)
    phats <- est_phats(Y, mrkrs, baseline, gamma, summary_fn = summary_fn)
    if (!is.null(references)) 
        phats <- phats[-unlist(pure_samples), ]
    
    return(list(estimates = phats, markers = mrkrs, n_markers = n_markers, gamma = gamma))
}
