#' Generate \code{assign} and \code{w.assign}.
#'
#' For every root node in DAG defined by \code{map}, \code{\link{paths}} circles
#' over all possible path graphs, picks up the one that consists of the
#' most unmarked node, and then marks the nodes in the path graph that have been
#' selected. \code{\link{paths}} won't move to the next root node, until all
#' the descendant nodes of the current root have been marked.
#'
#' @param map Matrix of \code{n.edges}-by-\code{2} dimension, where \code{n.edges}
#' is the number of directed edges in DAG. The first column has indices of
#' nodes that edges directing from, whereas the second column gives the indices
#' of nodes the corresponding edges directing towards. If a node indexed \code{i}
#' does not have edges linked to it, record the corresponding row as
#' \code{map[i, NA]}.
#' @param var Length-\code{n.nodes} list where \code{n.nodes} is the number
#' of nodes in DAG and for which the \code{l}th element contains
#' the indices of variables embedded in the \code{l}th node.
#' @param w Length-\code{n.nodes} vector of positive values for which
#' \code{w_l} gives the weight for \code{g_l}, where \code{n.nodes} is the number
#' of nodes in DAG. If this is \code{NULL}, \code{w_l = sqrt(|g_l|)} will be used.
#' Necessary condition is \code{w_l} increases with \code{|g_l|}.
#'
#' @return Returns \code{assign}, a matrix of \code{p} columns that gives the
#' assignments of variables over selected path graphs, and \code{w.assign}, a
#' list of the same length as the number of rows in \code{assign}.
#' \describe{
#' \item{\code{assign}: }{Each row of \code{assign} corresponds to a path graph decomposed from DAG.}
#' \item{\code{w.assign}: }{The \code{l}th element of the list contains the
#' weights corresponding to the \code{l}th row of \code{assign}
#' (the \code{l}th path graph).}
#' }
#'
#' @export
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
#' # Assume two parameters per node.
#' var <- as.list(data.frame(t(matrix(1:16, ncol=2, byrow=TRUE))))
#' paths.result <- paths(map, var)
#' paths.result$assign
#' paths.result$w.assign
paths <- function(map, var, w = NULL) {
  ## total number of nodes in the DAG
  n.nodes <- max(map, na.rm = TRUE)
  ## total number of variables embedded in nodes in the DAG
  p <- max(c(var, recursive = TRUE))
  stopifnot(length(var) == n.nodes)
  all.nodes <- 1:n.nodes
  roots <- all.nodes[!all.nodes %in% map[, 2]]
  selected.paths <- list()
  i <- 1  ## index for selected.paths
  ## record the nodes that have been marked
  marked <- NULL
  ## go through all the roots and record the chosen paths
  for (u in roots) {
    ## All the possible path graphs originated from node u
    all.paths <- path.find(u, map)
    ## all descendant nodes of node u
    all.descendants <- unique(c(all.paths, recursive = TRUE))
    while (!all(all.descendants %in% marked)) {
      if(!is.null(marked)) {
        all.paths <- sapply(all.paths,
          function(x) x[!x %in% marked], simplify = FALSE)
      }
      path.len <- sapply(all.paths, length)
      selected.paths[[i]] <- all.paths[[which.max(path.len)]]
      all.paths[[which.max(path.len)]] <- NULL
      marked <- c(marked, selected.paths[[i]])
      i <- i + 1
    }
  }
  ## generate the ancestor matrix
  ancestor <- matrix(0, ncol = n.nodes, nrow = n.nodes)
  for (u in 1:n.nodes) ancestor[u, ] <- ancestor.find(u, map, n.nodes)
  ## length of variables embedded in each node
  var.len <- sapply(var, length)
  ## generate the weight vector if it is NULL
  if (is.null(w)) {
    for (u in 1:n.nodes) w[u] <- sqrt(sum(var.len[which(ancestor[u, ] == 1)]))
  }
  ## generate the assign matrix and w.assign list
  n.paths <- length(selected.paths)
  assign <- matrix(0, nrow = n.paths, ncol = p)
  w.assign <- list()
  for (j in 1:n.paths) {
    w.assign[[j]] <- w[selected.paths[[j]]]
    last.var <- NULL
    for (k in 1:length(selected.paths[[j]])) {
      current.nodes <- which(ancestor[selected.paths[[j]][k], ] != 0)
      current.var <- rep(0, p)
      for (u in current.nodes) current.var[var[[u]]] <- 1
      if (is.null(last.var)) {
        indices <- current.var == 1
      } else {
        ## only record increment from last node
        indices <- current.var == 1 & last.var == 0
      }
      assign[j, indices] <- k
      last.var <- current.var
    }
  }
  return(list("assign" = assign, "w.assign" = w.assign))
}

#' Find ancestor nodes for a node in DAG.
#'
#' Recursively finds all ancestor nodes in DAG for node with the given index.
#'
#' @param index Index of the node currently at.
#' @param map Matrix of \code{n.edges}-by-\code{2} dimension, where \code{n.edges}
#' is the number of directed edges in DAG. The first column has indices of
#' nodes that edges directing from, whereas the second column gives the indices
#' of nodes the corresponding edges directing towards.
#' @param n.nodes Number of nodes in DAG.
#'
#' @return Returns a length-\code{n.nodes} vector of binary values for which
#' 1 indicates the corresponding node is an ancestor node and 0 indicates it is not.
ancestor.find <- function(index, map, n.nodes) {
  nodes <- rep(0, n.nodes)
  # Node with the given index is a root, whose ancestor node is itself.
  if (!index %in% map[, 2]) {
    nodes[index] <- 1
  # Trace up the DAG if the node with the given index is not a root.
  } else {
    left.index <- which(index == map[, 2])
    for (m in 1:length(left.index)) {
      nodes <- nodes + ancestor.find(map[left.index[m], 1], map, n.nodes)
    }
    nodes[c(index, which(nodes != 0))] <- 1
  }
  return(nodes)
}

#' Find all path graphs originated from a given root.
#'
#' Recursively find all possible path graphs originated from a given root in
#' DAG.
#'
#' @param index Index of a root node (a node whose index never appears in
#' \code{map[, 2]}).
#' @param map Matrix of \code{n.edges}-by-\code{2} dimension, where \code{n.edges}
#' is the number of directed edges in DAG. The first column has indices of
#' nodes that edges directing from, whereas the second column gives the indices
#' of nodes the corresponding edges directing towards.
#'
#' @return Returns a list of path graphs originated from root \code{index}, for
#' which the \code{i}th element of the returned list is a vector of indices of
#' nodes in the \code{i}th path graph.
path.find <- function(index, map) {
  # Returns index if the node is a leaf or the node is not linked to other nodes
  if (!index %in% map[, 1] || is.na(map[which(map[,1]==index), 2])) {
    return(list(index))
  } else {
    path <- list()
    ind <- 1
    right.index <- which(index == map[, 1])
    for (m in 1:length(right.index)) {
      path.right <- path.find(map[right.index[m], 2], map)
      for (n in 1:length(path.right)) {
        path[[ind]] <- c(index, path.right[[n]])
        ind <- ind + 1
      }
    }
    return(path)
  }
}
