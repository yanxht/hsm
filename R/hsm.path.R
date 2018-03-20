#' Solves proximal operator of latent group Lasso over a grid of lam values.
#'
#' See \code{\link{hsm}} for the problem that is solved. If \code{lamlist}
#' is not provided, a grid of lam values will be constructed starting at
#' \code{lammax}, the smallest value of lam for which the solution is completely
#' sparse.
#'
#' @param y Length-\code{p} vector used in proximal operator.
#' @param nlam Number of lam values to include in grid. Default value is 20.
#' @param flmin Ratio between the smallest lam and largest lam in grid. Default
#' value is 0.01. Increasing its value will give more sparse solutions.
#' @param lamlist A grid of lam values to use. If this is \code{NULL}, then
#' a grid of \code{nlam} lam values equally spaced in the logarithm scale
#' between \code{lammax} and \code{lammax * flmin} are used; otherwise,
#' \code{nlam} and \code{flmin} are ignored.
#' @param w Length-\code{n.nodes} vector of positive values for which
#' \code{w_l} gives the weight for \code{g_l}, where \code{n.nodes} is the number
#' of nodes in DAG. If this is \code{NULL},
#' \code{w_l = sqrt(|g_l|)} will be used.
#' Necessary condition is \code{w_l} increases with \code{|g_l|}.
#' @param map Matrix of \code{n.edges}-by-\code{2} dimension, where \code{n.edges}
#' is the number of directed edges in DAG. The first column has indices of
#' nodes that edges directing from, whereas the second column gives the indices
#' of nodes the corresponding edges directing towards. If a node indexed \code{i}
#' does not have edges linked to it, record the corresponding row as
#' \code{map[i, NA]}.
#' @param var Length-\code{n.nodes} list for which the \code{l}th element contains
#' the indices of variables embedded in the \code{l}th node.
#' @param assign Matrix of \code{p} columns that gives the assignments of
#' variables over different path graphs. Each row of \code{assign} corresponds to
#' a path graph decomposed from DAG. If this is NULL,
#' \code{\link{hsm}} first break down DAG into different path graphs,
#' and then give value to \code{assign} afterwards, based on \code{map} and
#' \code{var}; otherwise, \code{map} and \code{var} are ignored. Refer to
#' \code{\link{paths}} for more details.
#' @param w.assign List of length \code{nrow(assign)}, for which
#' the \code{l}th element contains the weights corresponding to the
#' \code{l}th row of \code{assign} (the \code{l}th path graph). For
#' example, if the \code{l}th path graph is made up of three nodes indexed with
#' \code{{3, 4, 6, 8}}, \code{w.assign[[l]] = {w_3, w_4, w_6, w_8}}. If this
#' is NULL, \code{\link{hsm}} will give value to \code{w.assign},
#' along with \code{assign}; otherwise, \code{map} and \code{var} are ignored.
#' Refer to \code{\link{paths}} for more details.
#' @param get.penalval If \code{TRUE}, \eqn{lam * \Omega(\beta; w)} are computed
#' and returned, otherwise \code{NA} is returned.
#' @param tol Tolerance level used in BCD. Convergence is assumed when no
#' parameter of interest in each path graph changes by more than tol in BCD.
#' @param maxiter Upperbound of the number of iterations that BCD to perform.
#'
#' @return Returns a sequence of estimates of the solution to the proximal
#' operator of the latent group Lasso. The returned solutions are exact ones if
#' the DAG is a directed path graph.
#' \item{lamlist}{Grid of lam values used.}
#' \item{beta.m}{A \code{nlam}-by-\code{p} matrix where
#' \code{beta.m[i, ]} gives the \code{i}th solution to the proximal
#' operator, corresponding to the \code{i}th lam value in the grid.}
#' \item{penalval.m}{Length-\code{nlam} vector of values of the penalty
#' \eqn{lam * \Omega(\beta; w)} where \code{penalval.m[i, ]} correponds to
#' the \code{i}th lam value in the grid, if \code{get.penalval} is \code{TRUE}.
#' If \code{get.penalval} is \code{FALSE}, \code{NA} is returned.}
#' \item{assign}{Value of \code{assign}.}
#' \item{w.assign}{Value of \code{w.assign}.}
#'
#' @export
#'
#' @seealso \code{\link{hsm}}
#' @seealso \code{\link{paths}}
#' @seealso \code{\link{lam.max.hsm}}
#'
#' @examples
#' # The following example appears in Figure 7 of Yan & Bien (2015).
#' # Generate map defining DAG.
#' map <- matrix(0, ncol=2, nrow=8)
#' map[1, ] <- c(1,2)
#' map[2, ] <- c(2,7)
#' map[3, ] <- c(3,4)
#' map[4, ] <- c(4,6)
#' map[5, ] <- c(6,7)
#' map[6, ] <- c(6,8)
#' map[7, ] <- c(3,5)
#' map[8, ] <- c(5,6)
#' # Assume one parameter per node.
#' # Let parameter and node share the same index.
#' var <- as.list(1:8)
#' set.seed(100)
#' y <- rnorm(8)
#' result <- hsm(y=y, lam=0.5, map=map, var=var, get.penalval=TRUE)
#' result.path <- hsm.path(y=y, map=map, var=var, get.penalval=TRUE)
hsm.path <- function(y, nlam = 20, flmin = 0.01,
  lamlist = NULL, w = NULL, map, var, assign = NULL, w.assign = NULL,
  get.penalval = FALSE, tol = 1e-8, maxiter = 1e4) {
  # Decompose the DAG into path graphs, if assign == NULL.
  if (is.null(assign)) {
    stopifnot(is.matrix(map))
    stopifnot(is.list(var))
    stopifnot(length(y) == max(c(var, recursive = TRUE), na.rm = TRUE))
    paths.result <- paths(map = map, var = var, w = w)
    assign <- paths.result$assign
    w.assign <- paths.result$w.assign
  }
  stopifnot(is.matrix(assign))
  stopifnot(is.list(w.assign))
  stopifnot(length(w.assign) == nrow(assign))
  stopifnot(length(y) == ncol(assign))
  ## Get a grid of lambda values
  if (is.null(lamlist)) {
    lammax <- lam.max.hsm(y, assign, w.assign)
    lamlist <- lammax * exp(seq(0, log(flmin), length = nlam))
  } else {
    nlam <- length(lamlist)
  }
  ## Compute over the grid of lambda values
  p <- length(y)
  beta.m <- matrix(0, nrow = nlam, ncol = p)
  penalval.m <- rep(0, length = nlam)
  for (i in seq(nlam)) {
    cat(i)
    result <- hsm(y = y, lam = lamlist[i], assign = assign,
      w.assign = w.assign, get.penalval = get.penalval, tol = tol,
      maxiter = maxiter)
    beta.m[i, ] <- result$beta
    if (get.penalval) penalval.m[i] <- result$penalval
  }
  if (!get.penalval) penalval.m <- NA
  return(list("lamlist" = lamlist, "beta.m" = beta.m,
    "penalval.m" = penalval.m, "assign" = assign, "w.assign" = w.assign))
}

#' Computes the smallest lam value such that beta = 0.
#'
#' Computes \code{lammax}, the smallest value of lam for which \code{\link{hsm}}
#' gives a completely sparse solution.
#'
#' @param y Length-\code{p} vector used in proximal operator.
#' @param assign Matrix of \code{p} columns that gives the assignments of
#' variables over different path graphs. Each row of \code{assign} corresponds to
#' a path graph decomposed from DAG. Refer to \code{\link{paths}} for more details.
#' @param w.assign List of length \code{nrow(assign)}, for which
#' the \code{l}th element contains the weights corresponding to the
#' \code{l}th row of \code{assign} (the \code{l}th path graph). For
#' example, if the \code{l}th path graph is made up of three nodes indexed with
#' \code{{3, 4, 6, 8}}, \code{w.assign[[l]] = {w_3, w_4, w_6, w_8}}.
#' Refer to \code{\link{paths}} for more details.

#'
#' @export
#'
#' @useDynLib hsm, .registration = TRUE
#'
#' @seealso \code{\link{hsm.path}}
lam.max.hsm <- function(y, assign, w.assign) {
  n.paths <- nrow(assign)
  p <- ncol(assign)
  lammaxs <- rep(0, length = n.paths)
  for (j in 1:n.paths) {
    out <- .C("pathgraph_lammax",
      as.double(y),
      as.double(w.assign[[j]]),
      as.integer(assign[j, ]),
      as.integer(p),
      as.integer(max(assign[j,])),
      lamval = as.double(lammaxs[j]))
    lammaxs[j] <- out$lamval
  }
 return(max(lammaxs))
}
