# ============================================
# = permutation test for the size of network =
# ============================================

network.test <- function(overall.network, gene.list, max.d = 1, n.sim = 1000) {
  obs <- grow.network(overall.network, gene.list, write.file = FALSE, max.distance = max.d)
  n.map <- length(intersect(V(obs$graph)$name, gene.list))
  n.max <- obs$max
  
  sim.max <- numeric(n.sim)
  
  for (i in 1:n.sim) {
    set.seed(i)
    sim.max[i] <- grow.network(overall.network, sample(V(overall.network)$name, n.map), write.file = FALSE, max.distance = max.d)$max
  }
  
  return( (sum(sim.max >= n.max) + 1)/(n.sim + 1) )
  
}