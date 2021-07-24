# ==============================================
# = function to perform r-spider like analysis =
# ==============================================

grow.network <- function(input.graph, seed, max.distance, file = "out.graphml", write.file = FALSE) {
  
  require(igraph)
  
  # make sure there are at least 2 genes in the seed list that intersect with the graph
  cat("input ", length(unique(seed)), " unique gene IDs.\n", sep = "")
  
  # calculate how many genes are mapped to the graph
  seed <- intersect(seed, V(input.graph)$name)
  cat(length(seed), " genes can be mapped to the global network.\n", sep = "")
  
  if (length(seed) < 2) {
    stop("fewer than 2 input nodes can be mapped to the global network")
  }
  
  # get neighbors of the seed, such that search for shortest paths can be smaller
  seed.neighbor.nodes <- unique(unlist(neighborhood(input.graph, order = max.distance - 1, nodes = V(input.graph)[V(input.graph)$name %in% seed], mode = "all")))
  input.graph <- simplify(induced.subgraph(input.graph, vids = seed.neighbor.nodes))
  
  # limit the search space to pairs of nodes whose shortest paths are within the threshold
  seed.node.index <- match(seed, V(input.graph)$name)
  short.path <- shortest.paths(input.graph, v = seed.node.index, to = seed.node.index, mode = "all")
  search.nodes <- which(short.path <= max.distance, arr.ind = T)
  
  # it is possible that none of the pairs are within the distance threshold
  output.graph <- graph.empty(directed = FALSE)
  if (nrow(search.nodes) == 0) {
    return(list(graph = output.graph, max = 0))
  }
  
  output.graph <- output.graph + V(input.graph)$name
  V(output.graph)$symbol <- V(input.graph)$symbol
  V(output.graph)$list <- "missing"
  V(output.graph)$list[seed.node.index] <- "list"
  
  for (i in 1:nrow(search.nodes)) {
    if (search.nodes[i, 1] != search.nodes[i, 2]) {
      add.paths <- get.all.shortest.paths(input.graph, from = seed.node.index[search.nodes[i, 1]], to = seed.node.index[search.nodes[i, 2]])$res
      for (j in 1:length(add.paths)) {
        output.graph <- output.graph + path(add.paths[[j]])
      }
    }
  }
  
  output.graph <- simplify(output.graph)
  output.degree <- degree(output.graph)
  output.graph <- induced.subgraph(output.graph, vids = unique(c(match(seed, V(output.graph)$name), which(output.degree >= 1))))

  # write to disk the graph in graphml format if instructed
  if (write.file) {
    write.graph(graph = output.graph, file = file, format = "graphml")
  }

  return(list(graph = output.graph, max = max(clusters(output.graph)$csize), mapped = length(seed)))
  
}