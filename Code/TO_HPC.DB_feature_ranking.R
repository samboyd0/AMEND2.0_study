#===============================#
#      Feature Ranking:          
# Degree Bias Mitigation Methods
#===============================#
# Clear environment
rm(list = ls())

# Path names
path.data = "/path/to/data/"
path.amend.code = "/path/to/AMEND/code/"
path.results = "/path/to/results/"

# Attach Libraries
library(igraph)
library(Matrix)
library(data.table)
library(biomaRt)

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
create_folds = function(kfolds, gene_sets){
  flag = unlist(lapply(gene_sets, length)) < kfolds
  if(any(flag)) stop(paste0("Not all pathways have at least kfolds=", kfolds, " features."))
  fold_ids = vector("list", length(gene_sets))
  names(fold_ids) = names(gene_sets)
  for(j in 1:length(gene_sets)){
    N = length(gene_sets[[j]])
    perm = sample(1:N, N)
    fold_ids[[j]] = cut(perm, breaks = min(kfolds,N), labels = FALSE) 
  }
  return(fold_ids)
}
getReactome <- function(species = 'human') {
  # Gets gene sets with Hugo gene symbols to use with FGSEA analysis.
  # species: currently accepts 'human' or 'mouse'
  require(reactome.db)
  require(annotate)
  db = ''
  if(species == 'human') {
    require(org.Hs.eg.db)
    db <- 'org.Hs.eg'
  } else if(species == 'mouse') {
    require(org.Mm.eg.db)
    db <- 'org.Mm.eg'
  } else {
    stop(paste0('Species ', species, ' not supported.'))
  }
  reactome_sets_full <- as.list(reactomePATHID2EXTID)
  pb <- txtProgressBar(min = 0, max = length(reactome_sets_full), style = 3)
  reactome_sets <- list()
  for(i in 1:length(reactome_sets_full)) {
    reactome_sets[[i]] <- as.vector(na.omit(getSYMBOL(reactome_sets_full[[i]], data=db)))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  names(reactome_sets) <- names(reactome_sets_full)
  xx = as.list(reactomePATHID2NAME)
  names(reactome_sets) = xx[names(reactome_sets)]
  reactome_sets <- reactome_sets[lapply(reactome_sets, length)>0]
  return(reactome_sets)
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

#================#
# Script Settings
#================#


path.S.drive = "/Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/Paper 3/Data/"
path.HPC = "/kuhpc/home/s703b506/AMEND2.0/DATA/"

path.ppi.el = "10090.protein.physical.links.v12.0.txt.gz"
path.ppi.g = "ppi_phys_graph.rds"
path.reactome = "mm_reactome_pathways.rds"
path.db.res = "DB_reactome_result_i"




#=====#
# DATA
#=====#
#=== Physical PPI Network ===#
edge_threshold = 0.7
ppi = data.table::fread(file = paste0(path.data, "10090.protein.physical.links.v12.0.txt.gz"), header = T) %>%
  dplyr::mutate(combined_score = combined_score / 1000) %>%
  dplyr::filter(combined_score >= edge_threshold) %>%
  dplyr::mutate(protein1 = extract_string(protein1, "\\.", 2), # remove taxonomy identifier
                protein2 = extract_string(protein2, "\\.", 2)) %>%
  dplyr::select(protein1, protein2, combined_score) %>%
  dplyr::distinct() %>%
  as.matrix()

# Remove duplicates and Aggregate edge scores
tmp = character(nrow(ppi))
for(i in seq_along(tmp)){
  tmp[i] = paste(sort(c(ppi[i,1], ppi[i,2])), collapse = "|")
}
scores = as.numeric(ppi[,3])
dat = stats::aggregate(scores, by = list(tmp), FUN = min)
ppi = matrix(c(extract_string(dat[,1], "\\|", 1), extract_string(dat[,1], "\\|", 2), dat[,2]), ncol = 3)

g = igraph::graph_from_edgelist(el = ppi[,1:2], directed = FALSE)
E(g)$weight = as.numeric(ppi[,3])

# Map ENSMUSP IDs to Gene Symbols
uniq.prot = unique(V(g)$name)
mart <- useMart(biomart = "ensembl", host = "https://jan2024.archive.ensembl.org")
mart_data <- useDataset("mmusculus_gene_ensembl", mart = mart)
# Attributes are what you want returned from your mapping (i.e., what you are mapping to)
# Filters are what you are trying to map 
# Example: If you want to map Symbols to Entrez_IDs, then filters=Symbols, attributes=c(Symbols, Entrez_IDs). attributes has both because we want to return both types.
prot.map <- getBM(attributes = c("ensembl_peptide_id", "mgi_symbol"),
                  filters = "ensembl_peptide_id",
                  values = uniq.prot,
                  mart = mart_data)
prot.map = prot.map[prot.map[,1] != "" & prot.map[,2] != "", ] # Remove features that didn't map
head(prot.map); nrow(prot.map)
# Check relationship of mapping 
table(table(prot.map[,1]))
table(table(prot.map[,2]))
# Remove genes that map to multiple proteins from prot.map
multi.gene = names(table(prot.map[,2]))[table(prot.map[,2]) > 1]
prot.map = prot.map[!prot.map[,2] %in% multi.gene,]
# Remove proteins from heterogeneous network that didn't map to a gene symbol
g = igraph::induced_subgraph(g, which(V(g)$name %in% prot.map[,1]))
# Map Symbols to nodes in network
id = match(V(g)$name, prot.map[,1])
igraph::vertex_attr(g, "mgi_symbol", which(!is.na(id))) = prot.map[id[!is.na(id)],2]
# Take largest connected component
g = largest_connected_component(g)

#=== Reactome Pathway Database ===#
if(0){
  r_paths = getReactome(species = "mouse")
  # Merging redundant pathways with "|" character
  tmp1 = unlist(lapply(r_paths, function(x) paste(sort(x), collapse = "|")))
  tmp2 = table(tmp1)
  tmp3 = vector("list", length(tmp2))
  for(i in seq_along(tmp3)){
    tmp3[[i]] = r_paths[[match(names(tmp2)[i], tmp1)]]
    nm = names(tmp1)[tmp1 %in% names(tmp2)[i]]
    nm = paste(nm, collapse = "|")
    names(tmp3)[i] = nm
    # names(tmp3)[i] = nm[which.min(nchar(nm))]
  }
  r_paths = tmp3
  # Removing genes not in PPIN
  g.symbols = V(g)$mgi_symbol
  r_paths = lapply(r_paths, function(x) x[x %in% g.symbols])
  rm.id = c()
  for(i in seq_along(r_paths)) if(length(r_paths[[i]]) == 0) rm.id = c(rm.id, i)
  r_paths = r_paths[-rm.id]
  # Save 
  saveRDS(r_paths, file = path.reactome)
}
r_paths = readRDS(paste0(path.data, "mm_reactome_pathways.rds"))

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
transition_matrix <- function(adjM, norm = c("degree", "core", "modified_degree"), k = 0.5){
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
# returns diffusion scores as a vector
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
# NB: Each method is deterministic.
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

#=====================#
# Feature Ranking with
# k-fold CV
#=====================#
# Save k-folds that were generated, or set seed for RNG
# For each transition matrix... (5)
#   For each adjustment method... (4)
#     For each pathway... (~1400) (<-- In Parallel)
#       Collect ranks of scores of test folds ('n.fold' runs)
#       Collect scores from using entire pathway as seeds (1 run)
# Total number of RWR runs: (n.fold + 1) * 1400 * 4 * 5 = 168,000 (n.fold=5)
# For figures: display the empirical CDF of ranks of nodes in hold-out folds. Plot ECDFs of different adj methods on same graph for a given transition matrix.

ll = 0.1
mm = 0.5
hh = 1

n.cores = 16

r.param = 0.5

adj.mat = igraph::as_adjacency_matrix(graph = g, attr = "weight", sparse = TRUE)
dimnames(adj.mat) = list(V(g)$mgi_symbol, V(g)$mgi_symbol)
node.names = rownames(adj.mat)

set.seed(40021514)
n.fold = 5
paths_cv = r_paths[unlist(lapply(r_paths, function(x) length(x) >= n.fold))] # Remove pathways with less than 'n.fold' features
fold.id = create_folds(kfolds = n.fold, gene_sets = paths_cv)

# t.mat.methods = c("degree", paste("modified_degree",ll,sep=";"), paste("modified_degree",mm,sep=";"), paste("modified_degree",hh,sep=";"), "core")
t.mat.methods = c("degree", paste("modified_degree",ll,sep=";"), paste("modified_degree",mm,sep=";"))
adj.methods = c("control", "tec", "bs", "in")
all.adj.fun = list(control, trans_ec_norm, bistochastic_scaling, inflate_normalize)
names(all.adj.fun) = adj.methods

# Degree Bias Analysis
tmat.window = list(1:3)[[1]]
adj.method.window = list(1:4)
l.part = 0
mid.point = 600

db.res = vector("list", length(t.mat.methods))
names(db.res) = t.mat.methods
for(i in tmat.window){ # seq_along(t.mat.methods)
  if(any(grepl("modified_degree", t.mat.methods[i]))){
    K = as.numeric(extract_string(t.mat.methods[i],";",2))
    t.mat.methods[i] = "modified_degree"
  }else K = NULL
  t.mat = transition_matrix(adjM = adj.mat, norm = t.mat.methods[i], k = K)
  db.res[[i]] = vector("list", length(adj.methods))
  names(db.res[[i]]) = adj.methods
  for(j in adj.method.window[[i]]){ # seq_along(adj.methods)
    adj.fun = all.adj.fun[[adj.methods[j]]]
    
    if(l.part == 1){
      parallel.iters = 1:mid.point
    }else if(l.part == 2){
      parallel.iters = (mid.point+1):length(paths_cv)
    }else if(l.part == 0){
      parallel.iters = 1:length(paths_cv)
    }
    
    cl = makeForkCluster(n.cores, outfile = "")
    registerDoParallel(cl)
    res = foreach(l = parallel.iters, .verbose = FALSE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dorng% {
      if(l %% 100 == 0) message(paste0("*** i",i,".j",j,".l",l," ***"))
      res.ii = matrix(0, ncol = n.fold+1, nrow = length(node.names))
      for(ii in seq_len(ncol(res.ii))){
        if(ii <= n.fold){ 
          train.set = paths_cv[[l]][!fold.id[[l]] %in% ii]
        }else train.set = paths_cv[[l]]
        seed.values = matrix(ifelse(node.names %in% train.set, 1, 0), ncol = 1, dimnames = list(node.names))
        res.ii[,ii] = adj.fun(trans_mat = t.mat, setSeeds = seed.values, restart = r.param)
      }
      res.ii
    }
    stopCluster(cl)
    db.res[[i]][[j]] = res
  }
}
# saveRDS(db.res, file = path.db.res)
# db.res = readRDS(path.db.res)

#=== Formating & Interpreting Results ===#
if(0){
  
  node.degree = Matrix::colSums(adj.mat)
  
  # db.res = vector("list", length(t.mat.methods)); names(db.res) = t.mat.methods
  for(i in seq_along(db.res)){
    # db.res[[i]] = vector("list", length(adj.methods)); names(db.res[[i]]) = adj.methods
    for(j in seq_along(db.res[[i]])){
      message(paste0("i",i,".j",j))
      # 'scores' will be same length & order as 'paths_cv'
      scores = readRDS(paste0(path.db.res, i, "j", j, ".rds"))
      # db.res[[i]][[j]] = vector("list", length(paths_cv))
      for(l in seq_along(db.res[[i]][[j]])){
        tmp.rank = c()
        tmp.degree = c()
        for(ii in seq_len(n.fold)){
          if(ii <= n.fold){
            test.set = paths_cv[[l]][fold.id[[l]] %in% ii]
          }else test.set = paths_cv[[l]]
          test.ids = which(node.names %in% test.set)
          # Look at ranks (lower ranks = better)
          if(ii <= n.fold){
            tmp.rank = c(tmp.rank, rank(1 / db.res[[i]][[j]][[l]][,ii])[test.ids])
          } 
        }
        db.res[[i]][[j]][[l]] = tmp.rank
      }
    }
  }
  
  #========#
  # Figures
  #========#
  # Save Dimensions
  # w/o legend: w=1687,h=1282
  # w/ legend: w=2500,h=1282
  # delta=30
  
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
  
  disc.color.theme = "Set1"
  div.color.theme = "RdBu"
  grad.color.theme ="Blues"
  

  id = 3
  t.mat.methods[id]
  ranks = unlist(db.res[[id]])
  ranks = c()
  groups = c()
  for(i in seq_along(db.res[[id]])){
    for(j in seq_along(db.res[[id]][[i]])){
      ranks = c(ranks, db.res[[id]][[i]][[j]])
      groups = c(groups, rep(adj.methods[i], length(db.res[[id]][[i]][[j]])))
    }
  }
  dat = data.frame(mr = ranks, Methods = groups, row.names = NULL)
  dat$Methods = factor(x = dat$Methods, levels = c('control', 'tec', 'in', 'bs'))
  levels(dat$Methods)
  # unique(sort(dat$Methods))
  plot_title = c("Degree Normalization", 
                 "Penalized Degree Normalization (k=0.1)",
                 "Penalized Degree Normalization (k=0.5)")[id]
  y_label = "Cumulative Probability"
  x_label = "Rank"
  delta = 30
  lw = 1.2
  
  # No legend
  if(0){
    p.db = ggplot(dat, aes(mr, colour = Methods)) +
      stat_ecdf(geom = "step", pad = FALSE, show.legend = TRUE, linewidth=lw) +
      labs(title = plot_title,
           x = x_label, y = y_label) +
      scale_color_manual(values = brewer.pal(n = 4, name = disc.color.theme)) + 
      # guides(color = "none") +
      theme(plot.title = element_text(size = 15+delta, face = "bold", hjust = 0),
            axis.title = element_text(size = 13+delta, face = "bold"),
            axis.text = element_text(size = 12+delta, face = "bold.italic"),
            legend.title = element_text(size = 12+delta, face = "bold.italic"),
            legend.text = element_text(size = 11+delta, face = "italic"),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major.y = element_line(colour = "black", linetype = "dashed", linewidth = 0.25))
    p.db
    if(1){
      x.lb = 0
      x.ub = c(200, 200)[id]
      y.ub = c(0.6, 0.6)[id]
      y.lb = 0
      inset = p.db +
        labs(title=NULL) +
        coord_cartesian(xlim = c(x.lb, x.ub), ylim = c(y.lb,y.ub)) +
        theme(legend.position="none",
              axis.title = element_text(size = 10+delta, face = "bold"),
              axis.text = element_text(size = 12+delta, face = "bold.italic"),
              panel.grid.major.y = element_blank())
      
      ggdraw(p.db + guides(color = "none")) +
        draw_plot(inset, x=0.35,y=0.15,width=0.5,height=0.5)
    }
  }
  # With Legend 
  if(0){
    p.db = ggplot(dat, aes(mr, colour = Methods)) +
      stat_ecdf(geom = "step", pad = FALSE, show.legend = TRUE, linewidth=lw) +
      labs(title = plot_title,
           x = x_label, y = y_label) +
      scale_color_manual(values = brewer.pal(n = 4, name = disc.color.theme), labels = c('Control', 'SDS', 'IN', 'BS')) + # labels = c("BS", "Control", "IN", "SDS") or c('Control', 'SDS', 'IN', 'BS')
      theme(plot.title = element_text(size = 15+delta, face = "bold", hjust = 0),
            axis.title = element_text(size = 13+delta, face = "bold"),
            axis.text = element_text(size = 12+delta, face = "bold.italic"),
            legend.title = element_text(size = 12+delta, face = "bold"),
            legend.text = element_text(size = 12+delta, face = "bold.italic"),
            legend.key.size = unit(1, units='cm'),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major.y = element_line(colour = "black", linetype = "dashed", linewidth = 0.25)) +
      guides(color = guide_legend(override.aes = list(linewidth = 4)))
    # p.db
    if(1){
      x.lb = 0
      x.ub = 200
      y.ub = 0.6
      y.lb = 0
      inset = p.db +
        labs(title=NULL) +
        coord_cartesian(xlim = c(x.lb, x.ub), ylim = c(y.lb,y.ub)) +
        theme(legend.position="none",
              axis.title = element_text(size = 10+delta, face = "bold"),
              axis.text = element_text(size = 12+delta, face = "bold.italic"),
              panel.grid.major.y = element_blank())
      
      ggdraw(p.db) +
        draw_plot(inset, x=0.2,y=0.15,width=0.5,height=0.5)
    }
  }
  
}


