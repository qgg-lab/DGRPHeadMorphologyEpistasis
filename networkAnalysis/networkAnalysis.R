# ==================================================
# = network analysis for a list of candidate genes =
# ==================================================

# 1. load genetic interaction data and necessary functions
# ============================================================

# specify what you want to do
# do you want missing genes?
n.missing.gene = 0

# libraries
library("igraph")

# functions
source("grow.network.R")
source("network.test.R")

# genetic interaction data in a RData
load("flybase.genet.int.RData")

# 2. now read the genes, prepare the data in two columns separated by tab
# first column is FBgn (required), second column is gene symbol, for plotting
# ============================================================

gene.list <- read.table("allGenes.txt", header = FALSE, as.is = TRUE)

# 3. extract network from the genetic interaction network
# max.distance = 1 means no missing gene
# max.distance = 2 means one missing gene, and so on
# ============================================================

genet.int.nomiss <- grow.network(flybase.genet.int.graph, gene.list[, 1], file = "allGenes.genet.nomiss.graphml", write.file = TRUE, max.distance = n.missing.gene + 1)
# not all genes would be on the global network
# the output graphml can be imported to cytoscape to draw networks
n.mapped.gene <- genet.int.nomiss$mapped

# 4. permutation test and calculate P value
# ============================================================

set.seed(1)
max.cluster <- numeric(1000)
for(i in 1:1000) {
  max.cluster[i] <- grow.network(flybase.genet.int.graph, sample(V(flybase.genet.int.graph)$name, n.mapped.gene), write.file = FALSE, max.distance = n.missing.gene + 1)$max
	cat(i, "\n")
}

# number of permutations with more than the max
sum(max.cluster >= genet.int.nomiss$max)
