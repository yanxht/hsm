#' Block coordinate descent based on path graphs for proximal operator of
#' latent group Lasso
#'
#' \code{hsm} is the R package implementing Algorithm 4 of Yan & Bien (2015)
#' that uses path-graph-based BCD to solve the proximal operator of latent
#' group Lasso in hierarchical sparse modeling (HSM). The algorithm solves
#' the proximal operator using BCD that circles over path graphs decomposed
#' a directed acyclic graph (DAG).
#'
#' The package is designed for situation in which latent group Lasso is used
#' to achieve hierarchical sparsity pattern in a DAG.
#' The hierarchical sparsity pattern is one when parameters embedded in a node being
#' set to zero, all the parameters embedded in the descendant nodes in DAG are
#' zeroed out as well.
#'
#' Its main functions are \code{\link{hsm}}, \code{\link{hsm.path}}.
#'
#' @author Xiaohan Yan \email{xy257@@cornell.edu}, Jacob Bien
#' @references Yan, X. and Bien, J. (2015) "Hierarchical
#' Sparse Modeling: A Choice of Two Regularizers"
#' @name hsm-package
#' @docType package
NULL
