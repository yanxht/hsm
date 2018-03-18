#' Solves proximal operator of latent group lasso in Hierarchical Sparse
#' Modeling.
#'
#' Solves proximal operator of the latent group Lasso appearing in Yan & Bien (2017)
#' \deqn{min_\beta || y - \beta ||_2^2 + lam * \Omega(\beta; w)}
#' where \eqn{\Omega(\beta; w) = min_{sum_l v^(l) = \beta; supp(v^(l))
#' \subset g_l} w_l * || v^(l) ||_2} is known as the latent group lasso penalty
#'  as defined in Jacob et al. (2009). In the problem, \eqn{\beta} is a length-\eqn{p}
#'  parameter vector and its elements are embedded in a directed acyclic graph
#'  (DAG). The desired sparsity pattern is a subgraph of the DAG such that if
#'  \eqn{\beta_i} embedded in node \eqn{i} are set to zero, all the parameters
#'  embedded in the descendant nodes of \eqn{i} are zeroed out as well. The problem
#'  is solved by breaking down the DAG into several path graphs for which closed-form
#'  solutions are available for the proximal operator corresponding with each
#'  path graph, and performing block coordinate descent across the path graphs.
#'  See Section 4.3 of the paper for more details and explanations.
#'
#' See Section 2.2 of the paper for problem setup and group structure specifications.
#' See Figure 7 in Section 4.3 for an example of decomposing DAG into path
#' graphs. See Algorithm 4 in paper for details of the path-based BCD.
#'
#' @param y Length-\code{p} vector used in proximal operator.
#' @param lam Non-negative tuning parameter that controls the sparsity level.
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
#' @param beta.ma \code{n.paths}-by-\code{p} matrix of initialization of
#' beta value in the \code{n.paths} path graphs decomposed from DAG.
#' Do not use unless you know the decomposition of DAG.
#'
#' @return Returns an estimate of the solution to the proximal operator of
#' the latent group Lasso. The returned value is an exact solution if the
#' DAG is a directed path graph.
#' \item{beta}{A length-\code{p} vector giving solution to the proximal
#' operator defined above.}
#' \item{ite}{Number of cycles of BCD performed.}
#' \item{penalval}{Value of the penalty \eqn{lam * \Omega(\beta; w)} if
#' \code{get.penalval} is \code{TRUE}, otherwise \code{NA}.}
#' \item{assign}{Value of \code{assign}.}
#' \item{w.assign}{Value of \code{w.assign}.}
#' \item{beta.ma}{\code{n.paths}-by-\code{p} matrix of decomposed beta
#' values for all the decomposed path graphs. The beta values are from the last iteration
#' in \code{\link{hsm}}.}
#'
#' @export
#'
#' @useDynLib hsm
#'
#' @references Yan, X. and Bien, J. (2017). Hierarchical Sparse Modeling:
#' A Choice of Two Group Lasso Formulations. Statist. Sci. 32,
#' no. 4, 531--560. doi:10.1214/17-STS622.
#' @references Jacob, L., Obozinski, G. and Vert, J. (2009). Group Lasso with Overlap and Graph Lasso.
#' In Proceedings of the 26th Annual International Conference on Machine Learning.
#' ICML'09 433-440. ACM, New York.
#'
#' @seealso \code{\link{hsm.path}}
#' @seealso \code{\link{paths}}
#'
#' @examples
#' # The following example appears in Figure 7 of Yan & Bien (2015).
#' # Generate map defining DAG.
#' map <- matrix(0, ncol=2, nrow=8)
#' map[1, ] <- c(1, 2)
#' map[2, ] <- c(2, 7)
#' map[3, ] <- c(3, 4)
#' map[4, ] <- c(4, 6)
#' map[5, ] <- c(6, 7)
#' map[6, ] <- c(6, 8)
#' map[7, ] <- c(3, 5)
#' map[8, ] <- c(5, 6)
#' # Assume one parameter per node.
#' # Let parameter and node share the same index.
#' var <- as.list(1:8)
#' set.seed(100)
#' y <- rnorm(8)
#' result <- hsm(y=y, lam=0.5, map=map, var=var, get.penalval=TRUE)
#'
#' # Another example in which DAG contains two separate nodes
#' map <- matrix(0, ncol=2, nrow=2)
#' map[1, ] <- c(1, NA)
#' map[2, ] <- c(2, NA)
#' # Assume ten parameters per node.
#' var <- list(1:10, 11:20)
#' set.seed(100)
#' y <- rnorm(20)
#' lam <- 0.5
#' result <- hsm(y=y, lam=lam, map=map, var=var, get.penalval=TRUE)
#' # The solution is equivalent to performing group-wise soft-thresholdings
#' beta.st <- c(y[1:10] * max(0, 1 - lam * sqrt(10) / norm(y[1:10], "2")),
#'           y[11:20] * max(0, 1 - lam * sqrt(10) / norm(y[11:20], "2")))
#' all.equal(result$beta, beta.st)
hsm <- function(y, lam, w = NULL, map, var, assign = NULL,
  w.assign = NULL, get.penalval = FALSE, tol = 1e-8, maxiter = 1e4,
  beta.ma = NULL) {
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
	# Do BCD on each path graph, represented by each row in assign.
  	p <- length(y)
  	# Number of decomposed path graphs.
	n.paths <- nrow(assign)
	beta0 <- beta1 <- NULL
	if (is.null(beta.ma)) {
	  beta0 <- beta1 <- matrix(0, ncol = p, nrow = n.paths)
	} else {
	  beta0 <- beta1 <- beta.ma
	}
	penalval <- rep(0, n.paths)
	continue <- TRUE
	ite <- 0
	while (continue) {
	  # If get.penalval == TRUE, compute the value of penalty.
	  if (get.penalval) {
	    for (i in 1:n.paths) {
	      out <- .C("pathgraph_prox2",
	        r = as.double(y - colSums(beta1[-i, ,drop = FALSE])),
	        as.double(lam),
	        as.double(w.assign[[i]]),
	        as.integer(assign[i,]),
	        as.integer(p),
	        as.integer(max(assign[i,])),
	        pen = as.double(penalval[i]))
	      beta1[i, ] <- out$r
	      penalval[i] <- out$pen
	    }
	  # No need to compute the value of penalty if get.penalval == FALSE.
	  } else {
	    for (i in 1:n.paths) {
	      out <- .C("pathgraph_prox",
	        r = as.double(y - colSums(beta1[-i, ,drop = FALSE])),
	        as.double(lam),
	        as.double(w.assign[[i]]),
	        as.integer(assign[i,]),
	        as.integer(p),
	        as.integer(max(assign[i,])))
	      beta1[i, ] <- out$r
	    }
	  }
		if (max(abs(beta0 - beta1)) < tol || ite > maxiter) continue <- FALSE
		beta0 <- beta1
		ite <- ite + 1
		#cat("Path-based BCD finishes ", ite, "iteration.", fill=TRUE)
	}
	return(list("beta" = colSums(beta1), "ite" = ite,
	  "penalval" = ifelse(get.penalval, sum(penalval), NA),
	  "assign" = assign, "w.assign" = w.assign, "beta.ma"=beta1))
}
