#' given a skid id and neuronlist return neuron name
#' @param nl a neuron list
#' @export
skid_name <- function (skid, nl) {
  nl[as.character(skid),'name']
}

#' Check occurrences of a tag in each neuron of a neuron list
#' @param nl a neuron list
#' @param tag name of the tag(e.g. 'pedunculus-calyx boundary')
#' @param counts number of occurrences of the tag in each neuron
#' @export
#' @return a list containing the number of tag occurrences in each neuron
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
#' @importFrom stats dist
soma_dist <- function (n1, n2) {
  dist(rbind(xyzmatrix(n1$d[n1$StartPoint,]), xyzmatrix(n2$d[n2$StartPoint,])))
}

#' get cable length of neuron
#' @details  return cable length given a skeleton id
#'    and a neuronlist that contains the neuron,
#' @export
cable_len <- function(nl, skid) {
  summary(nl[[as.character(skid)]])[['cable.length']]
}

#' a dfs search to traverse a neuron
#' @param starting_idx the node index of the starting point of the traversal
#' @export
walk_skeleton <- function(x, starting_idx) {
n_graph_dfs = igraph::dfs(as.ngraph(x), root = starting_idx, "out", unreachable = FALSE)$order
}

add_field <- function(n, field_content, field_name) {
  n[[field_name]] = field_content
  n
}

#' download neurons associated with an annotation
#' @details The annotation is attached in the meta-data in a field called anno_type.
#' Also include in each neuron a field called skids that contains skids. Skids serve
#' as a convenient id for the neuron.
#' @param anno an annotation used in CATMAID
#' @export
fetchn_by_annotation <- function(anno, conn, anno_type, ...) {
  nl = catmaid::read.neurons.catmaid(paste0("annotation:^", anno, "$"), conn=conn, ...)
  nl[, anno_type] = anno
  skid_ids = attr(nl, 'df')$skid
  nmapply(add_field, nl, skid_ids, 'skid')
}

#' download neurons associated with annotations
#' @description transform the neurons if ref is supplied,
#' process their anatomy segments if border_tag is supplied
#' @details currently only works if each neuron is annotated with only one of the annotations
#' @param annos a list of annotations used in CATMAID
#' @param border_tag border tag name (e.g 'bouton border' or 'claw border'),
#' if not (\code{NULL}) it will use the \code{\link{segment_arbor}} function to process the arbors
#' @export
#' @seealso
#' \code{\link{segment_arbor}}, \code{\link{fetchn_by_annotation}}
fetch_mb_neurons <- function(annos, conn, anno_type='annos', ref="FCWB", border_tag=NULL, ...) {
  nl = do.call(c, lapply(annos, fetchn_by_annotation, conn, anno_type, ...))

  if (anyDuplicated(nl[,"skid"]) != 0)  message("More than 1 annotations on a neuron!")
  if (!is.null(ref)) nl = nlapply(nl, xform_brain, sample="FAFB13", reference=ref)
  if (!is.null(border_tag)) nl = nlapply(nl, function(n) segment_arbor(n, border_tag, n$skid))
  nl
}


#' Find a neuron's anatomic segments delimited by border_tag
#' @description Anatomic segments could be boutons ('bouton border' tag)
#' or claws ('claw border' tag).
#' The tags label border nodes that define the anatomic segments in
#' alternate fashion.
#' Take boutons for example in dfs order,
#' Root -> 1st tagged node: non-bouton,
#' 1st tagged node -> 2nd tagged node: bouton,
#' 2nd tagged node -> 3rd tagged node: non-bouton
#' ...
#' @param n A neuron
#' @param border_tag name of the tag labeling the nodes that define the border
#' @param skid optional, but I always like to include skid as part of the unique
#' id for each anatomic segment.
#' @return neuron with 2 new fields: AnaSeg, a list with each element containing
#' node indices in a single anatmic segment (a bouton or a claw); BoutonEnds, indices of end nodes of
#' each terminal bouton
#' @export
segment_arbor <- function(n, border_tag, skid=NULL) {
  ng=as.ngraph(n)

  if (!"bouton border" %in% names(n$tags)) {
    n$AnaSeg = list()
    n$BoutonEnds = list()
    message(paste0("This neuron doesn't have tag ", border_tag))
  } else {

    # indices of all nodes with the tag
    idx=get_tagged_node(n, border_tag, "indices")

    # find incident edge to the tagged nodes, cut those edges,
    # find all connected components, each componnent would be a single
    # bouton/non-bouton segment
    border_vs=ego(ng, 1, nodes = idx, mode = "in")
    t5=do.call(c,lapply(border_vs, as_ids))
    c_border_vs=unlist(lapply(seq(t5), function(x) t5[[x + (x%%2)*2 -1]]))
    bb1_e=E(ng, c_border_vs)
    ng_d=delete_edges(ng, bb1_e)
    mshp=components(ng_d, mode = "weak")[["membership"]]

    # the component that contains the root (usually soma)
    seed_comp=mshp[[n$StartPoint]]

    rd_graph=make_graph(mshp[c_border_vs])

    # distance from the root to each AnaSeg
    path_to_root=shortest_paths(rd_graph, from = seed_comp, to = V(rd_graph), mode = "out")

    # boolean list of whether a component is bouton (claw) or not
    is_btn = sapply(path_to_root[[1]], function(x) 1-length(x)%%2)

    # assign all tagged nodes os part of anatomic segments (bouton or claw)
    for (i in idx) {
      comp = mshp[[i]]
      inci_id = setdiff(ego(rd_graph, 1, nodes = comp, mode = "in")[[1]], comp)
      if (is_btn[[comp]] == 0 & length(inci_id)==1) mshp[[i]]=inci_id
    }

    # ng = set_vertex_attr(ng, 'segment', value=mshp)
    # ng = set_vertex_attr(ng, 'is_bouton', value=is_btn[mshp])
    # seg_counts=sum(sapply(path_to_root[[1]], function(x) 1-length(x)%%2))

    # add bouton endpoints
    btn_ends = c()
    for (i in n$EndPoints) {
      if (is_btn[[mshp[[i]]]]==1) btn_ends=c(btn_ends,i)
    }

    dist_to_root = sapply(path_to_root[[1]], function(x) length(x)-1)

    # gives an id to each anatomic segment
    # the id consisits of skid (or _). group_boolean . ids in reduced_graph . distance to root
    ana_seg = list()
    for (i in seq_along(dist_to_root)) {
      if (is_btn[[i]]==1) {
        btn_id = paste(ifelse(is.null(skid),"_",skid), is_btn[[i]], i, dist_to_root[[i]], sep = '.')
        ana_seg[[btn_id]] = which(mshp==i)
      }
    }
    n$AnaSeg = ana_seg
    n$BoutonEnds = btn_ends
  }
  n
}
