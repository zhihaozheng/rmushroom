#' given a skid id and neuronlist return neuron name
#' @param nl a neuron list
#' @export
skid_name <- function (skid, nl) {
  attr(nl[as.character(skid)],'df')$name
}

#' check occurances of a tag in each neuron of a neuron list
#' @param nl a neuron list
#' @param tag name of the tag(e.g. 'pedunculus-calyx boundary')
#' @param counts number of occurances of the tag in each neuron
#' @export
#' @return a list containing the number of tag occurances in each neuron
check_tag <- function(nl, tag, counts=1) {
  c=nlapply(nl, function(nl) length(nl$tags[[tag]]))
  if (sum(c==counts) != length(c)) warning ('Tag counts per neuron do not match!')
  c
}

#' return attributes of nodes with a given tag
#'
#' @param n a neuron object
#' @param tag name of the tag(e.g. 'pedunculus-calyx boundary')
#' @param return_fields type of output to return
#'   ("indices" or any combination of c("PointNo","Label","X","Y","Z","W","Parent"))
#' @export
get_tagged_node <- function(n, tag, rval) {
  nodes = n$tags[[tag]]
  if (length(nodes)==0) stop('No such tag!')
  if (length(rval)==1 & rval=='indices') return(which(n$d$PointNo %in% nodes))
  n$d[n$d$PointNo %in% nodes, rval]
}

#' return an unique node index around a section number
#' @details given a neuron and a section number, return a unique node
#'    index in the neuron closest to the section. Note it will look
#'    up and down until it finds an unique node (exclude no node or
#'    more than 1 node)
#' @export
z_to_node_idx <- function (n, z) {
  z = z*35
  idx = which(n$d$Z == z)
  m = 2
  while (length(idx) != 1) {
    z = z + (1 - (m%%2)*2)*35*(m %/% 2)
    m = m + 1
    idx = which(n$d$Z == z)
  }
  return (idx)
}

#' return nearest neighbour distance for a pair of neurons
#' @export
get_nndist <- function(query, target) {
  nnres = nabor::knn(xyzmatrix(target), xyzmatrix(query), k = 1)
  return(drop(nnres$nn.dists))
}

#' return soma distance for a pair of neurons
#' @export
soma_dist <- function (n1, n2) {
  dist(rbind(xyzmatrix(n1$d[n1$StartPoint,]), xyzmatrix(n2$d[n2$StartPoint,])))
}

#' get cable length of neuron
#' @details  return cable length given a skeleton id
#'    and a neuronlist that contains the neuron,
#' @export
cable_len <- function(nl, skid) {
  unlist(summary(nl[as.character(skid)], include.attached.dataframe = FALSE)['cable.length'])
}

#' a dfs search to traverse a neuron
#' @param starting_idx the node index of the starting point of the traversal
#' @export
walk_skeleton <- function(x, starting_idx) {
n_graph_dfs = igraph::dfs(as.ngraph(x), root = starting_idx, "out", unreachable = FALSE)$order
}
