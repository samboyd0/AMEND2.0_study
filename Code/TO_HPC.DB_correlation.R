#===============================#
#      Degree Correlation:          
# Degree Bias Mitigation Methods
#===============================#
# Clear environment
rm(list = ls())

# Path names
path.data = "/path/to/data/"
path.amend.code = "/path/to/AMEND/code/"
path.results = "/path/to/results/"
path.cor.res = paste0(path.results, "DB_correlation_results.rds")

# Attach Libraries
library(igraph)
library(Matrix)
library(data.table)

library(doParallel)
library(foreach)
library(doRNG)

# Helper functions
extract_string = function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))
largest_connected_component = function(g){
  if(!igraph::is_connected(g)){
    message("Taking largest connected component.")
    comps = igraph::components(g)
    largest_comp_id = which.max(comps$csize)
    g = igraph::induced_subgraph(g, which(comps$membership == largest_comp_id))
  }
  g
}
get_diagonal = function(X){
  if(!"dgCMatrix" %in% class(X)) X = Matrix::Matrix(X, sparse = TRUE)
  d = numeric(ncol(X))
  for(i in 1:ncol(X)){
    id.tmp = (X@p[i] + 1):X@p[i+1] # ids of X@x that are non-zero and in col i
    row.ids = X@i[id.tmp] + 1 # row ids of non-zero elements in col i
    if(i %in% row.ids){ # if diagonal is non-zero
      d[i] = X@x[id.tmp[row.ids == i]]
    }else next
  }
  d
}

#=====#
# DATA
#=====#
# Load in DE analysis results (from get_Omics.R)
omics = list(als = readRDS(paste0(path.data, "ALS_DE_results.rds")),
             lc = readRDS(paste0(path.data, "LC_DE_results.rds")),
             uc = readRDS(paste0(path.data, "UC_DE_results.rds")),
             cd = readRDS(paste0(path.data, "CD_DE_results.rds")),
             hd = readRDS(paste0(path.data, "HD_DE_results.rds")))
omics = lapply(omics, function(x){
  x$entrez_id = as.numeric(rownames(x))
  x
})

# Load in PPIN (from get_string_igraph.R)
ppi = readRDS(paste0(path.data, "STRING_0.5_igraph.rds"))

# Remove genes from GE data that don't map to PPIN
# lapply(omics, function(x) sum(x$entrez_id %in% V(ppi)$name))
omics = lapply(omics, function(x) x[x$entrez_id %in% V(ppi)$name,])

# Add data as vertex attributes to ppi
for(i in seq_along(omics)){
  ids = match(V(ppi)$name, omics[[i]]$entrez_id)
  igraph::vertex_attr(ppi, paste(names(omics)[i], "lfc", sep="_"), which(!is.na(ids))) = omics[[i]]$logFC[ids[!is.na(ids)]]
  igraph::vertex_attr(ppi, paste(names(omics)[i], "pval", sep="_"), which(!is.na(ids))) = omics[[i]]$P.Value[ids[!is.na(ids)]]
}

#========#
# METHODS
#========#
sum2one <- function(X) {
  x.dimnames = dimnames(X)
  if(!"dgCMatrix" %in% class(X)) X = methods::as(methods::as(methods::as(X, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  inv_col_sum = Matrix::Diagonal(x = Matrix::colSums(X)^(-1))
  res = X %*% inv_col_sum
  dimnames(res) = x.dimnames
  res[is.na(res)] = 0 # This ensures that columns with all zeros that were divided by zero remain zero.
  res
}
stationary.distr = function(x){
  e = Re(RSpectra::eigs(A = x, k = 1, which = "LM")$vectors[,1])
  tmp = e / sum(e)
  ifelse(tmp < 0, 0, tmp)
}
entropy = function(x){
  tmp = ifelse(round(x, 10) == 0, 0, -x * log(x))
  sum(tmp)
}

#=== Transition Matrix ===#
# 3 options: degree, biased degree, coreness
transition_matrix <- function(adjM, norm = c("degree", "core", "modified_degree", "maximum_entropy"), k = 0.5){
  # Coerce to a general, column-compressed sparse matrix from Matrix package
  if(!"dgCMatrix" %in% class(adjM)) adjM = methods::as(methods::as(methods::as(adjM, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  norm <- match.arg(norm)
  if(norm == "core"){
    g = igraph::graph_from_adjacency_matrix(adjmatrix = adjM, mode = "undirected", weighted = TRUE)
    core = igraph::coreness(g)
    nadjM = adjM
    for(j in 1:ncol(nadjM)){
      id.tmp = adjM@i[(adjM@p[j] + 1):adjM@p[j+1]] + 1 # row ids of non-zero elements in column j
      nadjM[id.tmp, j] = core[id.tmp] / sum(core[id.tmp])
    }
  }else if(norm == "degree"){
    wt <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-1))
    nadjM <- adjM %*% wt
  }else if(norm == "modified_degree"){
    wt <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-k))
    adjM.mod <- wt %*% adjM %*% wt
    wt.mod <- Matrix::Diagonal(x = Matrix::colSums(adjM.mod)^(-1))
    nadjM <- adjM.mod %*% wt.mod
  }else{
    nadjM <- adjM
  }
  dimnames(nadjM) = dimnames(adjM)
  return(nadjM)
}

#=== RWR ===#
RWR <- function (nadjM, setSeeds = NULL, restart = 0.5){
  if (is.null(restart) || is.na(restart) || restart < 0 || restart > 100) {
    r <- 0.75
  }else if (restart > 1 && restart < 100) {
    r <- restart/100
  }else {
    r <- restart
  }
  stop_delta <- 1e-06
  stop_step <- 50
  if (is.null(setSeeds)) {
    P0matrix <- Matrix::Matrix(diag(nrow(nadjM)), sparse = T)
    rownames(P0matrix) <- rownames(nadjM)
    colnames(P0matrix) <- rownames(nadjM)
  }else {
    if (is.matrix(setSeeds) | is.data.frame(setSeeds)) {
      data <- as.matrix(setSeeds)
    }else if (is.vector(setSeeds)) {
      data <- matrix(setSeeds, ncol = 1, dimnames = list(names(setSeeds)))
    }
    if (is.null(rownames(data))) {
      stop("The function must require the row names of the input setSeeds.\n")
    }else if (any(is.na(rownames(data)))) {
      warning("setSeeds with NA as row names will be removed")
      data <- data[!is.na(rownames(data)), ]
    }
    ind <- match(rownames(data), rownames(nadjM)) # This ensures that all(rownames(data) == rownames(nadjM)[ind]) is TRUE
    nodes_mapped <- rownames(nadjM)[ind[!is.na(ind)]]
    if (length(nodes_mapped) != nrow(nadjM)) {
      warning("The row names of input setSeeds do not contain all those in the input graph.\n")
    }
    P0matrix <- matrix(0, nrow = nrow(nadjM), ncol = ncol(data))
    P0matrix[ind[!is.na(ind)], ] <- as.matrix(data[!is.na(ind),]) # ensuring that order of P0matrix is same as rownames(nadjM)
    P0matrix <- sum2one(P0matrix)
    # P0matrix <- Matrix::Matrix(P0matrix, sparse = T)
  }
  if (restart == 1) {
    PTmatrix <- P0matrix
  }else {
    PTmatrix <- Matrix::Matrix(0, nrow = nrow(P0matrix), ncol = ncol(P0matrix), sparse = T)
    P0 <- P0matrix[, 1]
    step <- 0
    delta <- 1
    PT <- P0
    while (delta > stop_delta && step <= stop_step) {
      PX <- (1 - r) * nadjM %*% PT + r * P0
      delta <- sum(abs(PX - PT))
      PT <- PX
      step <- step + 1
    }
    PT[PT < 1e-10] <- 0
    PTmatrix[, 1] <- PT
  }
  PTmatrix <- sum2one(PTmatrix)
  PTmatrix[PTmatrix < 1e-10] <- 0
  rownames(PTmatrix) <- rownames(P0matrix)
  colnames(PTmatrix) <- colnames(P0matrix)
  res = as.vector(PTmatrix)
  names(res) = rownames(PTmatrix)
  res
}

#==========================#
#=== Adjustment Methods ===#
#==========================#
# 4 options: control, Bistochastic scaling, Inflation-Normalization, Scaling by Stationary Distribution (Eigenvector centrality of transition matrix)
# Each will take as input an adjacency matrix or a transition matrix and seed values, and it will return a vector of importance scores
#=== Bistochastic Scaling ===#
# Iterative Proportional Fitting
ipf = function(X, e = 1e-6, gamma = 1){
  # For starting transition matrix X, obtain B = PXQ, where B is bistochastic, P,Q are diagonal matrices, and B and X are as similar as possible (minimize relative entropy between B & X)
  # Set target row sums depending on how much to homogenize degree influence. gamma = 1 for complete homogenization, gamma = 0 for no correction.
  get_rowsum_targets = function(rs, l, t) rs - l * (rs - t)
  target.row.sums = get_rowsum_targets(Matrix::rowSums(X), gamma, 1)
  target.col.sums = 1
  
  x.dimnames = dimnames(X)
  
  # Adding self-loops to aid in convergence
  d = get_diagonal(X)
  d[d == 0] = e
  X = X + Matrix::Diagonal(n = nrow(X), x = d)
  
  stop_delta = 1e-6
  step = 1
  stop_step = 200
  q1 = rep(1, nrow(X)) # initialize the diagonal elements of Q
  p1 = rep(1, nrow(X)) # initialize the diagonal elements of P
  repeat{
    # set p, given q (p are the diagonal elements of P)
    p2 = target.row.sums / Matrix::rowSums(X %*% Matrix::Diagonal(x = q1))
    # set q, given p found above (q are the diagonal elements of Q)
    q2 = target.col.sums / Matrix::colSums(Matrix::Diagonal(x = p2) %*% X)
    
    delta.p = all(abs(p2 - p1) <= stop_delta)
    delta.q = all(abs(q2 - q1) <= stop_delta)
    step = step + 1
    if((delta.p && delta.q) || step > stop_step) break
    q1 = q2
    p1 = p2
  }
  P = Matrix::Diagonal(x = p2)
  Q = Matrix::Diagonal(x = q2)
  B = P %*% X %*% Q
  dimnames(B) = x.dimnames
  return(list(B = B, p = p2, q = q2))
}
bistochastic_scaling = function(trans_mat, gamma = 1, setSeeds = NULL, restart = 0.5, ...){
  if(gamma == 0) return(trans_mat)
  x = 1e-6
  if(!"dgCMatrix" %in% class(trans_mat)) trans_mat = methods::as(methods::as(methods::as(trans_mat, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  B.tmp = ipf(trans_mat, x, gamma)
  B = B.tmp$B
  b = get_diagonal(B)
  tmp.res = numeric(Matrix::nnzero(B))
  for(i in 1:ncol(B)){
    if(b[i] == 0) next
    id.tmp = (B@p[i] + 1):B@p[i+1] # ids of B@x that are non-zero and in col i
    row.ids = B@i[id.tmp] + 1 # row ids of non-zero elements in col i
    diag.id = id.tmp[row.ids == i]
    off.id = id.tmp[row.ids != i]
    # Evenly distribute to neighbors of node i. Preserves column sums
    degr = B@p[i+1] - B@p[i] - 1 # number of non-zero elements in col i i.e., degree of node i. Minus 1 b/c of self-loops
    tmp = b[i] / degr
    tmp.res[diag.id] = 0
    tmp.res[off.id] = B@x[off.id] + tmp
  }
  B@x = tmp.res
  B = Matrix::drop0(B)
  # return(B)
  # Run RWR
  scores = RWR(B, setSeeds, restart)
  return(scores)
}

#=== Inflation-normalization ===#
inflate_normalize <- function (trans_mat, setSeeds = NULL, restart = 0.5, ...){
  # Changing to a Row-compressed sparse matrix
  # trans_mat = methods::as(as.matrix(trans_mat), "dgRMatrix")
  trans_mat = methods::as(methods::as(methods::as(trans_mat, "dMatrix"), "generalMatrix"), "RsparseMatrix")
  # Getting stationary distribution of trans_mat
  stat.distr1 = stationary.distr(trans_mat)
  e0 = entropy(stat.distr1)
  # Performing Inflation/normalization on transition matrix
  kf = c(1, 10, seq(50, 2000, 50))
  res = numeric(length(kf))
  for(j in seq_along(kf)){
    inflation = 1 + kf[j] * stat.distr1
    trans_mat.tmp = trans_mat
    tmp = numeric(Matrix::nnzero(trans_mat.tmp))
    for(i in seq_along(inflation)){ # i corresponds to rows of trans_mat
      # id.tmp = (trans_mat.tmp@p[i] + 1):trans_mat.tmp@p[i+1] # ids of trans_mat.tmp@x that are non-zero and in row i
      # n.tmp = trans_mat.tmp@p[i+1] - trans_mat.tmp@p[i] # number of non-zero elements in row i
      # tmp[id.tmp] = rep(inflation[i], n.tmp)
      tmp[(trans_mat.tmp@p[i] + 1):trans_mat.tmp@p[i+1]] = rep(inflation[i], trans_mat.tmp@p[i+1] - trans_mat.tmp@p[i])
    }
    trans_mat.tmp@x = trans_mat.tmp@x ^ tmp
    trans_mat.tmp = sum2one(trans_mat.tmp)
    sd.tmp = stationary.distr(trans_mat.tmp)
    if(any(sd.tmp < 0 | sd.tmp > 1)){
      if(j == 1){
        trans_mat = methods::as(methods::as(methods::as(trans_mat, "dMatrix"), "generalMatrix"), "CsparseMatrix")
        # return(trans_mat)
        scores = RWR(trans_mat, setSeeds, restart)
        return(scores)
      }else{
        j = j-1
        break
      }
    }
    res[j] = entropy(sd.tmp)
    if(j == 1){
      if(res[j] < e0){
        trans_mat = methods::as(methods::as(methods::as(trans_mat, "dMatrix"), "generalMatrix"), "CsparseMatrix")
        # return(trans_mat)
        scores = RWR(trans_mat, setSeeds, restart)
        return(scores)
      }
    }else if(res[j] <= res[j-1]){
      break
    }
  }
  res = res[1:j]
  j = which.max(res)
  inflation = 1 + kf[j] * stat.distr1
  tmp = numeric(Matrix::nnzero(trans_mat))
  for(i in seq_along(inflation)){ # i corresponds to rows of trans_mat
    tmp[(trans_mat@p[i] + 1):trans_mat@p[i+1]] = rep(inflation[i], trans_mat@p[i+1] - trans_mat@p[i])
  }
  trans_mat@x = trans_mat@x ^ tmp
  trans_mat = sum2one(trans_mat)
  # return(trans_mat)
  # Run RWR
  scores = RWR(trans_mat, setSeeds, restart)
  return(scores)
}

#=== Stationary Distribution Scaling ===#
# Otherwise known as eigenvector centrality scaling (eigenvector centrality of transition matrix, not adjacency matrix)
trans_ec_norm = function(trans_mat, setSeeds = NULL, restart = 0.5, ...){
  raw_scores = RWR(trans_mat, setSeeds, restart)
  ec = stationary.distr(trans_mat)
  # ec = RWR(nadjM = trans_mat, restart = 0)
  ec = ifelse(ec == 0, min(ec[ec != 0]), ec) # To avoid dividing by zero
  adjusted_scores = raw_scores / ec
  return(adjusted_scores / sum(adjusted_scores))
}

#=== Control ===#
control = function(trans_mat, setSeeds = NULL, restart = 0.5, ...){
  raw_scores = RWR(trans_mat, setSeeds, restart)
  return(raw_scores)
}

#============================#
# Degree Correlation Analysis ----
#============================#
# For each transition matrix... (5)
#   For each adjustment method... (4)
#     For each GE dataset... (5)
#       Collect scores using p-values as seeds (1 run)
#       Collect scores using |logFC| as seeds (1 run)
# Total number of RWR runs: 2 * 5 * 4 * 5 = 160


ll = 0.1
mm = 0.5
hh = 1

n.cores = 4

r.param = 0.5

adj.mat = igraph::as_adjacency_matrix(graph = ppi, attr = "weight", sparse = TRUE)
node.names = rownames(adj.mat)

t.mat.methods = c("degree", paste("modified_degree",ll,sep=";"), paste("modified_degree",mm,sep=";"), paste("modified_degree",hh,sep=";"), "core")
adj.methods = c("control", "tec", "bs", "in")
all.adj.fun = list(control, trans_ec_norm, bistochastic_scaling, inflate_normalize)
names(all.adj.fun) = adj.methods

# Degree Bias Correlation Analysis
# NB: No adj method or transition matrix type will change the order of the matrix.
db.res = vector("list", length(t.mat.methods))
names(db.res) = t.mat.methods
for(i in seq_along(t.mat.methods)){
  message(paste0("i",i))
  if(any(grepl("modified_degree", t.mat.methods[i]))){
    K = as.numeric(extract_string(t.mat.methods[i],";",2))
    t.mat.methods[i] = "modified_degree"
  }else K = NULL
  t.mat = transition_matrix(adjM = adj.mat, norm = t.mat.methods[i], k = K)
  db.res[[i]] = vector("list", length(adj.methods))
  names(db.res[[i]]) = adj.methods
  for(j in seq_along(adj.methods)){
    message(paste0("j",j))
    adj.fun = all.adj.fun[[adj.methods[j]]]
    cl = makeForkCluster(n.cores, outfile = "")
    registerDoParallel(cl)
    res = foreach(l = 1:length(omics), .verbose = FALSE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dorng% {
      res.ii = matrix(0, ncol = 2, nrow = length(node.names), dimnames = list(node.names, c("pval", "lfc")))
      for(ii in seq_len(ncol(res.ii))){
        v.nm = paste(names(omics)[l],colnames(res.ii)[ii],sep="_")
        seed.values = vertex_attr(ppi, v.nm)
        if(colnames(res.ii)[ii] == "pval") seed.values = -log10(seed.values + 1e-6) else seed.values = abs(seed.values)
        seed.values[is.na(seed.values)] = 0
        seed.values = matrix(seed.values, ncol = 1, dimnames = list(node.names))
        res.ii[,ii] = adj.fun(trans_mat = t.mat, setSeeds = seed.values, restart = r.param)
      }
      res.ii
    }
    stopCluster(cl)
    db.res[[i]][[j]] = res
  }
}
# saveRDS(db.res, file = path.cor.res)
# db.res = readRDS(path.cor.res)

#========#
# Figures ----
#========#
if(0){
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
  
  disc.color.theme = "Set1"
  div.color.theme = "RdBu"
  grad.color.theme ="Blues"
  
  t.mat.ids = 1:3
  
  node.degree = Matrix::colSums(adj.mat)
  node.cor.deg = vector("list", length(t.mat.ids))
  names(node.cor.deg) = names(db.res)[t.mat.ids]
  for(i in t.mat.ids){
    node.cor.deg[[i]] = vector("list", length(db.res[[i]]))
    names(node.cor.deg[[i]]) = names(db.res[[i]])
    for(j in seq_along(adj.methods)){
      node.cor.deg[[i]][[j]] = vector("list", length(db.res[[i]][[j]]))
      names(node.cor.deg[[i]][[j]]) = names(db.res[[i]][[j]])
      for(l in seq_along(omics)){
        node.cor.deg[[i]][[j]][[l]] = numeric(2)
        for(ii in seq_len(2)){
          node.cor.deg[[i]][[j]][[l]][ii] = cor(node.degree, db.res[[i]][[j]][[l]][,ii])
        }
        names(node.cor.deg[[i]][[j]][[l]]) = colnames(db.res[[i]][[j]][[l]])
      }
    }
  }
  
  #=== Figures ===#
  # https://bookdown.org/skaltman/visualization-book/discrete-continuous.html#two-discrete-categories
  # https://stulp.gmw.rug.nl/ggplotworkshop/twodiscretevariables.html
  
  d.t = c("pval", "lfc")[2] # seed 
  sum.fun = mean
  val = group = c()
  for(i in seq_along(node.cor.deg)){ # tmat
    for(j in seq_along(node.cor.deg[[i]])){ # adj method
      tmp = unlist(lapply(node.cor.deg[[i]][[j]], function(x) x[d.t]))
      val = c(val, sum.fun(tmp)) # take summary of correlations
      tmp1 = ifelse(grepl("modified",names(node.cor.deg)[i]), paste0(gsub(";"," (",names(node.cor.deg)[i]),")"), names(node.cor.deg)[i])
      group = c(group, paste(tmp1,names(node.cor.deg[[i]])[j],sep="|"))
    }
  }
  # sort(unique(extract_string(group, "\\|", 2)))
  dat = data.frame(val = val,
                   group = group,
                   tmat = factor(extract_string(group, "\\|", 1), levels=c('degree', 'modified_degree (0.1)', 'modified_degree (0.5)')),
                   adjmeth = factor(extract_string(group, "\\|", 2), levels=c("control", "tec", 'in', 'bs')))
  
  col.tmp = brewer.pal(n = length(unique(dat$adjmeth)), name = disc.color.theme)

  delta = 30
  
  ggplot(dat, aes(x = tmat, y = val, fill = adjmeth)) +
    geom_col(position = "dodge", width=0.7) +
    scale_x_discrete(labels = c("Degree", "Penalized Degree (0.1)", "Penalized Degree (0.5)")) +
    scale_fill_manual(values = col.tmp, labels = c('Control', 'SDS', 'IN', 'BS')) +
    labs(title = "Correlation of Diffusion Scores with Node Degree",
         subtitle = paste0("Seed Type: ", ifelse(d.t == 'pval', '-log10(P-Value)', '|LogFC|')),
         x = "Transition Matrix",
         y = "Pearson Correlation",
         fill = "Adjustment Method") +
    coord_cartesian(ylim=c(-0.25,0.7)) + 
    theme(plot.title = element_text(size=15+delta,hjust=0,face='bold'),
          plot.subtitle = element_text(size=10+delta,hjust=0,face='bold'),
          axis.title = element_text(size=15+delta,face="bold"),
          legend.title = element_text(size=12+delta,face="bold"),
          legend.text = element_text(size=10+delta,face='bold.italic'),
          axis.text.x.bottom = element_text(size=10+delta,face='bold.italic',angle=15,vjust=0.52,margin=margin(t=5,b=10)),
          axis.text.y.left = element_text(size=14+delta,face='bold.italic',margin=margin(r=10)),
          axis.ticks.x = element_line(linewidth=5),
          axis.ticks.length.x = unit(x=0.5,units='cm'),
          panel.background = element_blank(),
          panel.border = element_rect(color='lightgrey',fill=NA,linetype=1,linewidth=0.7),
          panel.grid.minor.y = element_line(color='lightgrey',linetype=2,linewidth=0.2),
          panel.grid.major.y = element_line(color='lightgrey',linetype=1,linewidth=0.7))
    
  
  # Save Dimensions
  # w=1800, h=1344
  # delta=25
  
  
}
