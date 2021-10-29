library("igraph")
library("readr")
library("colorspace")

## DCM vs HCM
network <- read.delim(file = "dcm_vs_hcm/network_dcm_vs_hcm_100.txt")
attributes <- read.delim(file = "dcm_vs_hcm/attributes_dcm_vs_hcm_100.txt")

df <- as.data.frame(network[, c(1, 3)])
gg <- graph_from_data_frame(d = df, directed = TRUE)
adj <- get.adjacency(graph = gg)

# Add node attributes
myrhombus <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color,
          stars = cbind(1.2*vertex.size, vertex.size, 1.2*vertex.size, vertex.size),
          add = TRUE, inches = FALSE)
}

add_shape("rhombus", clip = shapes("circle")$clip,
          plot = myrhombus)

ascii <- diverge_hcl(n = 201)

nodes <- V(gg)
col_type <- c()
node_type <- c()
for(ii in 1:length(nodes)){
  
  currNode <- rownames(adj)[nodes[ii]]
  idx <- which(attributes[, 1]==currNode)
  
  if(attributes[idx, 6]==""){
    node_type <- c(node_type, "circle")
  } else {
    if(attributes[idx, 6]=="T"){
      node_type <- c(node_type, "square")
    } else {
      node_type <- c(node_type, "rhombus")
    }
  }
  
  if(currNode=="Perturbation"){
    col_type <- c(col_type, "#00FF00")
  } else {
    col_type <- c(col_type, ascii[as.numeric(attributes[idx, 5])+101])
  }
  
}

V(gg)$color <- col_type
V(gg)$shape <- node_type
plot.igraph(gg)
dev.off()

# Add edge attributes
edges <- E(gg)
col_type <- c()
for(ii in 1:nrow(network)){
  
  if(network[ii, 2]=="1"){
    col_type <- c(col_type, "#000000")
  } else {
    col_type <- c(col_type, "#CC0000")
  }
}

E(gg)$color <- col_type

pdf(file = "dcm_vs_hcm/igraph_plot_dcm_vs_hcm.pdf", width = 10, height = 10)
plot.igraph(gg)
dev.off()

save(gg, file = "dcm_vs_hcm/igraph_dcm_vs_hcm.RData")

## DCM vs NFD
network <- read.delim(file = "dcm_vs_nfd/network_dcm_vs_nfd_100.txt")
attributes <- read.delim(file = "dcm_vs_nfd/attributes_dcm_vs_nfd_100.txt")

df <- as.data.frame(network[, c(1, 3)])
gg <- graph_from_data_frame(d = df, directed = TRUE)
adj <- get.adjacency(graph = gg)

# Add node attributes
myrhombus <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color,
          stars = cbind(1.2*vertex.size, vertex.size, 1.2*vertex.size, vertex.size),
          add = TRUE, inches = FALSE)
}

add_shape("rhombus", clip = shapes("circle")$clip,
          plot = myrhombus)

ascii <- diverge_hcl(n = 201)

nodes <- V(gg)
col_type <- c()
node_type <- c()
for(ii in 1:length(nodes)){
  
  currNode <- rownames(adj)[nodes[ii]]
  idx <- which(attributes[, 1]==currNode)
  
  if(attributes[idx, 6]==""){
    node_type <- c(node_type, "circle")
  } else {
    if(attributes[idx, 6]=="T"){
      node_type <- c(node_type, "square")
    } else {
      node_type <- c(node_type, "rhombus")
    }
  }
  
  if(currNode=="Perturbation"){
    col_type <- c(col_type, "#00FF00")
  } else {
    col_type <- c(col_type, ascii[as.numeric(attributes[idx, 5])+101])
  }
  
}

V(gg)$color <- col_type
V(gg)$shape <- node_type
plot.igraph(gg)
dev.off()

# Add edge attributes
edges <- E(gg)
col_type <- c()
for(ii in 1:nrow(network)){
  
  if(network[ii, 2]=="1"){
    col_type <- c(col_type, "#000000")
  } else {
    col_type <- c(col_type, "#CC0000")
  }
}

E(gg)$color <- col_type

pdf(file = "dcm_vs_nfd/igraph_plot_dcm_vs_nfd.pdf", width = 10, height = 10)
plot.igraph(gg)
dev.off()

save(gg, file = "dcm_vs_nfd/igraph_dcm_vs_nfd.RData")

## HCM vs NFD
network <- read.delim(file = "hcm_vs_nfd/network_hcm_vs_nfd_100.txt")
attributes <- read.delim(file = "hcm_vs_nfd/attributes_hcm_vs_nfd_100.txt")

df <- as.data.frame(network[, c(1, 3)])
gg <- graph_from_data_frame(d = df, directed = TRUE)
adj <- get.adjacency(graph = gg)

# Add node attributes
myrhombus <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color,
          stars = cbind(1.2*vertex.size, vertex.size, 1.2*vertex.size, vertex.size),
          add = TRUE, inches = FALSE)
}

add_shape("rhombus", clip = shapes("circle")$clip,
          plot = myrhombus)

ascii <- diverge_hcl(n = 201)

nodes <- V(gg)
col_type <- c()
node_type <- c()
for(ii in 1:length(nodes)){
  
  currNode <- rownames(adj)[nodes[ii]]
  idx <- which(attributes[, 1]==currNode)
  
  if(attributes[idx, 6]==""){
    node_type <- c(node_type, "circle")
  } else {
    if(attributes[idx, 6]=="T"){
      node_type <- c(node_type, "square")
    } else {
      node_type <- c(node_type, "rhombus")
    }
  }
  
  if(currNode=="Perturbation"){
    col_type <- c(col_type, "#00FF00")
  } else {
    col_type <- c(col_type, ascii[as.numeric(attributes[idx, 5])+101])
  }
  
}

V(gg)$color <- col_type
V(gg)$shape <- node_type
plot.igraph(gg)
dev.off()

# Add edge attributes
edges <- E(gg)
col_type <- c()
for(ii in 1:nrow(network)){
  
  if(network[ii, 2]=="1"){
    col_type <- c(col_type, "#000000")
  } else {
    col_type <- c(col_type, "#CC0000")
  }
}

E(gg)$color <- col_type

pdf(file = "hcm_vs_nfd/igraph_plot_hcm_vs_nfd.pdf", width = 10, height = 10)
plot.igraph(gg)
dev.off()

save(gg, file = "hcm_vs_nfd/igraph_hcm_vs_nfd.RData")
