#' Blocklist of genes for dimensionality reduction
#' 
#' A list of signatures, for mouse and human. These include cell cycling,
#' heat-shock genes, mitochondrial genes, and other genes classes, that may
#' confound the identification of cell types. These are used internally by scGate
#' and excluded from the calculation of dimensional reductions (PCA).
#' 
#' @name genes.blacklist.default
#' @docType data
#' @format A list of signatures
NULL

#' Toy dataset to test the package
#' 
#' A downsampled version (300 cells) of the single-cell dataset by Zilionis et
#' al. (2019) <doi:10.1016/j.immuni.2019.03.009>, with precalculated PCA and
#' UMAP reductions.
#' 
#' @name query.seurat
#' @docType data
#' @format A Seurat object
NULL
