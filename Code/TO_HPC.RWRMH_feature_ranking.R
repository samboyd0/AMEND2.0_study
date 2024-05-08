#============#
### TO HPC ###
#============#
rm(list=ls())

# Attach libraries
library(igraph)
library(Matrix)
library(doParallel)
library(foreach)
library(doRNG)

# Helper Functions
create_folds = function(kfolds, gene_sets){
  fold_ids = vector("list", length(gene_sets))
  names(fold_ids) = names(gene_sets)
  for(j in 1:length(gene_sets)){
    N = length(gene_sets[[j]])
    perm = sample(1:N, N)
    fold_ids[[j]] = cut(perm, breaks = min(kfolds,N), labels = FALSE) 
  }
  return(fold_ids)
}

#============#
# LOAD DATA
#============#
TO_HPC = FALSE
ppi_edge_threshold = 0.9
mmi_edge_threshold = 0.9

if(TO_HPC){
  g.het1 = readRDS(paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/het_net1.rds"))
  g.het2 = readRDS(paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/het_net2.rds"))
  g.mh = readRDS(paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/mh_net.rds"))
  g.multi.homo = readRDS(paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/multi.homo_net.rds"))
  
  kp.h1 = readRDS(paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10,"/kp_h1.rds"))
  kp.h2 = readRDS(paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10,"/kp_h2.rds"))
  kp.mh = readRDS(paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10,"/kp_mh.rds"))
  kp.multi.homo = readRDS(paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10,"/kp_multi.homo.rds"))

  source("/kuhpc/home/s703b506/AMEND2.0/CODE/amend_code/run_AMEND.R")
  source("/kuhpc/home/s703b506/AMEND2.0/CODE/amend_code/create_integrated_graph.R")
  source("/kuhpc/home/s703b506/AMEND2.0/CODE/amend_code/utils.R")
  source("/kuhpc/home/s703b506/AMEND2.0/CODE/amend_code/RandomWalk.R")
}else{
  g.het1 = readRDS(paste0("/Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/Paper 3/Data/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/het_net1.rds"))
  g.het2 = readRDS(paste0("/Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/Paper 3/Data/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/het_net2.rds"))
  g.mh = readRDS(paste0("/Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/Paper 3/Data/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/mh_net.rds"))
  g.multi.homo = readRDS(paste0("/Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/Paper 3/Data/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/multi.homo_net.rds"))
  
  kp.h1 = readRDS(paste0("/Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/Paper 3/Data/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10,"/kp_h1.rds"))
  kp.h2 = readRDS(paste0("/Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/Paper 3/Data/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10,"/kp_h2.rds"))
  kp.mh = readRDS(paste0("/Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/Paper 3/Data/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10,"/kp_mh.rds"))
  kp.multi.homo = readRDS(paste0("/Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/Paper 3/Data/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10,"/kp_multi.homo.rds"))
  
  source("/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Code/amend_code/run_AMEND.R")
  source("/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Code/amend_code/RandomWalk.R")
  source("/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Code/amend_code/create_integrated_graph.R")
  source("/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Code/amend_code/utils.R")
}

#===============================#
# Feature Ranking with k-fold CV ----
#===============================#
# For each transition matrix... (4)
#   For each Kegg Pathway... (~1000) (<-- In Parallel)
#     For each fold... (k.folds)
#       Get ranks of nodes in hold-out fold (1)

k.folds = 5
n.cores = 12
# Parameters for RWR
ll = 0.15
mm = 0.5
hh = 0.85

#===============================#
#=== Multiplex-heterogeneous ===#
#===============================#
tmat.res.path = paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/transition_matrix_multiplex_heterogeneous.rds")
rank.res.path = paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/node.ranking.results_multiplex_heterogeneous.rds")
if(0){
  message("*** Multi-Hetero ***")
  
  # Creating adjacency matrices
  adj.mat.mh = igraph::as_adjacency_matrix(g.mh, attr = "weight", sparse = TRUE)
  
  # Parameters for RWR
  tmat.params = list(low = list(jp = c(meta = ll, prot = ll),
                                sl = list(prot = c(prot_phys = ll, prot_fun = ll))),
                     medium = list(jp = c(meta = mm, prot = mm),
                                   sl = list(prot = c(prot_phys = mm, prot_fun = mm))),
                     high = list(jp = c(meta = hh, prot = hh),
                                 sl = list(prot = c(prot_phys = hh, prot_fun = hh))),
                     mixed = list(jp = c(meta = hh, prot = hh),
                                  sl = list(prot = c(prot_phys = ll, prot_fun = ll))))
  seed.params = list(nw = c(meta = 0.5, prot = 0.5),
                     lw = list(prot = c(prot_phys = 0.5, prot_fun = 0.5)))
  # Construct transition matrices
  if(1){
    tmat.mh = vector("list", 5); names(tmat.mh) = c("low", "medium", "high", 'mixed', "not_mh")
    for(i in seq_along(tmat.mh)){
      message(i)
      if(i == length(tmat.mh)){
        tmat.mh[[i]] = transition_matrix(adjM = adj.mat.mh, norm = "degree", heterogeneous = FALSE, multiplex = FALSE)
      }else{
        tmat.mh[[i]] = transition_matrix(adjM = adj.mat.mh, norm = "degree", heterogeneous = TRUE, multiplex = TRUE,
                                                jump.prob = tmat.params[[i]]$jp, switch.layer.prob = tmat.params[[i]]$sl)
      }
    }
    saveRDS(tmat.mh, file = tmat.res.path)
  }else tmat.mh = readRDS(tmat.res.path)
  
  # Aggregation function
  agg.fun = mean
  # Restart parameter
  r.param = 0.5
  
  #=== k-fold CV ===#
  set.seed(565)
  kegg_cv = kp.mh[unlist(lapply(kp.mh, function(x) length(x$GENE) >= k.folds && length(x$COMPOUND) >= k.folds))]
  gene_cv = lapply(kegg_cv, function(x) unname(x$GENE))
  comp_cv = lapply(kegg_cv, function(x) unname(x$COMPOUND))
  gene_fold_ids = create_folds(kfolds = k.folds, gene_sets = gene_cv)
  comp_fold_ids = create_folds(kfolds = k.folds, gene_sets = comp_cv)
  
  # Not Parallel
  if(0){
    tmat.res = vector("list", length(tmat.mh)); names(tmat.res) = names(tmat.mh)
    for(i in seq_along(tmat.mh)){ # For each transition matrix (low, high, mono-homo)
      message(paste0("i",i))
      node.names = extract_string(rownames(tmat.mh[[i]]), "\\|", 1)
      tmat.res[[i]] = vector("list", length(kegg_cv))
      for(j in seq_along(kegg_cv)){ # For each pathway
        message(paste0("j",j))
        tmat.res[[i]][[j]] = vector("list", k.folds)
        for(l in seq_len(k.folds)){ # For each fold
          # message(paste0("l",l))
          # Access features of fold l for set j 
          train.set = c(gene_cv[[j]][!gene_fold_ids[[j]] %in% l],
                        comp_cv[[j]][!comp_fold_ids[[j]] %in% l])
          test.set = c(gene_cv[[j]][gene_fold_ids[[j]] %in% l],
                       comp_cv[[j]][comp_fold_ids[[j]] %in% l])
          # Create matrix of seed values
          Seeds = matrix(ifelse(node.names %in% train.set, 1, 0), ncol = 1, dimnames = list(rownames(tmat.mh[[i]])))
          if(names(tmat.mh)[i] == "not_mh") mh.flag = FALSE else mh.flag = TRUE
          PTmatrix2 = RWR(nadjM = tmat.mh[[i]], setSeeds = Seeds, restart = r.param, heterogeneous = mh.flag, multiplex = mh.flag, net.weight = seed.params$nw, layer.weight = seed.params$lw)
          rwr_score = as.vector(PTmatrix2)
          names(rwr_score) = orig.names
          # Aggregate RWR scores
          tmp.val = stats::aggregate(rwr_score, by = list(agg.names), FUN = agg.fun)
          rwr_score = tmp.val$x
          names(rwr_score) = extract_string(tmp.val$Group.1, "\\|", 1)
          # Get median rank of test fold
          ranks = rank(rwr_score)
          tmat.res[[i]][[j]][[l]] = ranks[names(rwr_score) %in% test.set]
        }
      }
    }
  }
  # Parallel
  if(1){
    tmat.res = vector("list", length(tmat.mh)); names(tmat.res) = names(tmat.mh)
    for(i in seq_along(tmat.mh)){ # For each transition matrix (low, mid, high, mixed, mono-homo)
      message(paste0("i",i))
      orig.names = rownames(tmat.mh[[i]])
      agg.names = V(g.mh)$agg_name[match(orig.names, V(g.mh)$name)]
      node.names = extract_string(orig.names, "\\|", 1)
      cl = makeForkCluster(n.cores, outfile = "")
      registerDoParallel(cl)
      res = foreach(j = 1:length(kegg_cv), .verbose = TRUE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dorng% {
        tmp.res = vector("list", k.folds)
        for(l in seq_len(k.folds)){ # For each fold
          # message(paste0("l",l))
          # Access features of fold l for set j 
          train.set = c(gene_cv[[j]][!gene_fold_ids[[j]] %in% l],
                        comp_cv[[j]][!comp_fold_ids[[j]] %in% l])
          test.set = c(gene_cv[[j]][gene_fold_ids[[j]] %in% l],
                       comp_cv[[j]][comp_fold_ids[[j]] %in% l])
          # Create matrix of seed values
          Seeds = matrix(ifelse(node.names %in% train.set, 1, 0), ncol = 1, dimnames = list(orig.names))
          if(names(tmat.mh)[i] == "not_mh") mh.flag = FALSE else mh.flag = TRUE
          PTmatrix2 = RWR(nadjM = tmat.mh[[i]], setSeeds = Seeds, restart = r.param, heterogeneous = mh.flag, multiplex = mh.flag, net.weight = seed.params$nw, layer.weight = seed.params$lw)
          rwr_score = as.vector(PTmatrix2)
          names(rwr_score) = orig.names
          # Aggregate RWR scores
          tmp.val = stats::aggregate(rwr_score, by = list(agg.names), FUN = agg.fun)
          rwr_score = tmp.val$x
          names(rwr_score) = extract_string(tmp.val$Group.1, "\\|", 1)
          # Get ranks of test fold
          ranks = rank(1 / rwr_score)
          tmp.res[[l]] = ranks[names(rwr_score) %in% test.set]
        }
        tmp.res
      }
      stopCluster(cl)
      tmat.res[[i]] = res
    }
  }
  saveRDS(tmat.res, file = rank.res.path)
}

#===============================#
#=== Multiplex-Homogeneous ===#
#===============================#
tmat.res.path = paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/transition_matrix_multiplex_homogeneous.rds")
rank.res.path = paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/node.ranking.results_multiplex_homogeneous.rds")
if(0){
  message("*** Multi-Homo ***")
  
  # Creating adjacency matrices
  adj.mat.mh = igraph::as_adjacency_matrix(g.multi.homo, attr = "weight", sparse = TRUE)
  
  # Parameters for RWR
  tmat.params = list(low = list(sl = list(prot = c(prot_phys = ll, prot_fun = ll))),
                     medium = list(sl = list(prot = c(prot_phys = mm, prot_fun = mm))),
                     high = list(sl = list(prot = c(prot_phys = hh, prot_fun = hh))))
  seed.params = list(lw = list(prot = c(prot_phys = 0.5, prot_fun = 0.5)))
  # Construct transition matrices
  if(0){
    tmat.mh = vector("list", 5); names(tmat.mh) = c("low", "medium", "high", 'mixed', "not_mh")
    for(i in seq_along(tmat.mh)){
      message(i)
      if(i == length(tmat.mh)){
        tmat.mh[[i]] = transition_matrix(adjM = adj.mat.mh, norm = "degree", heterogeneous = FALSE, multiplex = FALSE)
      }else{
        tmat.mh[[i]] = transition_matrix(adjM = adj.mat.mh, norm = "degree", heterogeneous = FALSE, multiplex = TRUE,
                                         switch.layer.prob = tmat.params[[i]]$sl)
      }
    }
    saveRDS(tmat.mh, file = tmat.res.path)
  }else tmat.mh = readRDS(tmat.res.path)
  
  # Aggregation function
  agg.fun = mean
  # Restart parameter
  r.param = 0.5
  
  #=== k-fold CV ===#
  set.seed(565)
  kegg_cv = kp.multi.homo[unlist(lapply(kp.multi.homo, function(x) length(x$GENE) >= k.folds))]
  gene_cv = lapply(kegg_cv, function(x) unname(x$GENE))
  # comp_cv = lapply(kegg_cv, function(x) unname(x$COMPOUND))
  gene_fold_ids = create_folds(kfolds = k.folds, gene_sets = gene_cv)
  # comp_fold_ids = create_folds(kfolds = k.folds, gene_sets = comp_cv)
  
  # Not Parallel
  if(0){
    tmat.res = vector("list", length(tmat.mh)); names(tmat.res) = names(tmat.mh)
    for(i in seq_along(tmat.mh)){ # For each transition matrix (low, high, mono-homo)
      message(paste0("i",i))
      node.names = extract_string(rownames(tmat.mh[[i]]), "\\|", 1)
      tmat.res[[i]] = vector("list", length(kegg_cv))
      for(j in seq_along(kegg_cv)){ # For each pathway
        message(paste0("j",j))
        tmat.res[[i]][[j]] = vector("list", k.folds)
        for(l in seq_len(k.folds)){ # For each fold
          # message(paste0("l",l))
          # Access features of fold l for set j 
          train.set = c(gene_cv[[j]][!gene_fold_ids[[j]] %in% l],
                        comp_cv[[j]][!comp_fold_ids[[j]] %in% l])
          test.set = c(gene_cv[[j]][gene_fold_ids[[j]] %in% l],
                       comp_cv[[j]][comp_fold_ids[[j]] %in% l])
          # Create matrix of seed values
          Seeds = matrix(ifelse(node.names %in% train.set, 1, 0), ncol = 1, dimnames = list(rownames(tmat.mh[[i]])))
          if(names(tmat.mh)[i] == "not_mh") mh.flag = FALSE else mh.flag = TRUE
          PTmatrix2 = RWR(nadjM = tmat.mh[[i]], setSeeds = Seeds, restart = r.param, heterogeneous = mh.flag, multiplex = mh.flag, net.weight = seed.params$nw, layer.weight = seed.params$lw)
          rwr_score = as.vector(PTmatrix2)
          names(rwr_score) = orig.names
          # Aggregate RWR scores
          tmp.val = stats::aggregate(rwr_score, by = list(agg.names), FUN = agg.fun)
          rwr_score = tmp.val$x
          names(rwr_score) = extract_string(tmp.val$Group.1, "\\|", 1)
          # Get median rank of test fold
          ranks = rank(rwr_score)
          tmat.res[[i]][[j]][[l]] = ranks[names(rwr_score) %in% test.set]
        }
      }
    }
  }
  # Parallel
  if(1){
    tmat.res = vector("list", length(tmat.mh)); names(tmat.res) = names(tmat.mh)
    for(i in seq_along(tmat.mh)){ # For each transition matrix (low, mid, high, mono-homo)
      message(paste0("i",i))
      orig.names = rownames(tmat.mh[[i]])
      agg.names = V(g.multi.homo)$agg_name[match(orig.names, V(g.multi.homo)$name)]
      node.names = extract_string(orig.names, "\\|", 1)
      cl = makeForkCluster(n.cores, outfile = "")
      registerDoParallel(cl)
      res = foreach(j = 1:length(kegg_cv), .verbose = TRUE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dorng% {
        tmp.res = vector("list", k.folds)
        for(l in seq_len(k.folds)){ # For each fold
          # message(paste0("l",l))
          # Access features of fold l for set j 
          train.set = c(gene_cv[[j]][!gene_fold_ids[[j]] %in% l])
          test.set = c(gene_cv[[j]][gene_fold_ids[[j]] %in% l])
          # Create matrix of seed values
          Seeds = matrix(ifelse(node.names %in% train.set, 1, 0), ncol = 1, dimnames = list(orig.names))
          if(names(tmat.mh)[i] == "not_mh") mh.flag = FALSE else mh.flag = TRUE
          PTmatrix2 = RWR(nadjM = tmat.mh[[i]], setSeeds = Seeds, restart = r.param, heterogeneous = FALSE, multiplex = mh.flag, layer.weight = seed.params$lw)
          rwr_score = as.vector(PTmatrix2)
          names(rwr_score) = orig.names
          # Aggregate RWR scores
          tmp.val = stats::aggregate(rwr_score, by = list(agg.names), FUN = agg.fun)
          rwr_score = tmp.val$x
          names(rwr_score) = extract_string(tmp.val$Group.1, "\\|", 1)
          # Get ranks of test fold
          ranks = rank(1 / rwr_score)
          tmp.res[[l]] = ranks[names(rwr_score) %in% test.set]
        }
        tmp.res
      }
      stopCluster(cl)
      tmat.res[[i]] = res
    }
  }
  saveRDS(tmat.res, file = rank.res.path)
}

#================================#
#=== Monoplex-heterogeneous 1 ===#
#================================#
tmat.res.path = paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/transition_matrix_monoplex_heterogeneous1.rds")
rank.res.path = paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/node.ranking.results_monoplex_heterogeneous1.rds")
if(0){
  message("*** H1 ***")
  # Creating adjacency matrices
  adj.mat.h = igraph::as_adjacency_matrix(g.het1, attr = "weight", sparse = TRUE)
  dimnames(adj.mat.h) = list(paste(rownames(adj.mat.h), V(g.het1)$node_type, sep = "|"),
                             paste(rownames(adj.mat.h), V(g.het1)$node_type, sep = "|"))
  
  # Parameters for RWR
  tmat.params = list(low = list(jp = c(meta = ll, prot = ll)),
                     medium = list(jp = c(meta = mm, prot = mm)),
                     high = list(jp = c(meta = hh, prot = hh)))
  seed.params = list(nw = c(meta = 0.5, prot = 0.5))
  # Construct transition matrices
  if(1){
    tmat.h = vector("list", 4); names(tmat.h) = c("low", "medium", "high", "not_mh")
    for(i in seq_along(tmat.h)){
      message(i)
      if(i == length(tmat.h)){
        tmat.h[[i]] = transition_matrix(adjM = adj.mat.h, norm = "degree", heterogeneous = FALSE, multiplex = FALSE)
      }else{
        tmat.h[[i]] = transition_matrix(adjM = adj.mat.h, norm = "degree", heterogeneous = TRUE, multiplex = FALSE,
                                        jump.prob = tmat.params[[i]]$jp)
      }
    }
    saveRDS(tmat.h, file = tmat.res.path)
  }else tmat.h = readRDS(tmat.res.path)
  
  # Aggregation function
  agg.fun = mean
  # Restart parameter
  r.param = 0.5
  
  #=== k-fold CV ===#
  set.seed(565)
  kegg_cv = kp.h1[unlist(lapply(kp.h1, function(x) length(x$GENE) >= k.folds && length(x$COMPOUND) >= k.folds))]
  gene_cv = lapply(kegg_cv, function(x) unname(x$GENE))
  comp_cv = lapply(kegg_cv, function(x) unname(x$COMPOUND))
  gene_fold_ids = create_folds(kfolds = k.folds, gene_sets = gene_cv)
  comp_fold_ids = create_folds(kfolds = k.folds, gene_sets = comp_cv)
  
  tmat.res = vector("list", length(tmat.h)); names(tmat.res) = names(tmat.h)
  for(i in seq_along(tmat.h)){ # For each transition matrix (low, high, mono-homo)
    message(paste0("i",i))
    node.names = extract_string(rownames(tmat.h[[i]]), "\\|", 1)
    # Not parallel
    if(0){
      tmat.res[[i]] = vector("list", length(kegg_cv))
      for(j in seq_along(kegg_cv)){ # For each pathway
        message(paste0("j",j))
        tmat.res[[i]][[j]] = vector("list", k.folds)
        for(l in seq_len(k.folds)){ # For each fold
          # message(paste0("l",l))
          # Access features of fold l for set j 
          train.set = c(gene_cv[[j]][!gene_fold_ids[[j]] %in% l],
                        comp_cv[[j]][!comp_fold_ids[[j]] %in% l])
          test.set = c(gene_cv[[j]][gene_fold_ids[[j]] %in% l],
                       comp_cv[[j]][comp_fold_ids[[j]] %in% l])
          # Create matrix of seed values
          Seeds = matrix(ifelse(node.names %in% train.set, 1, 0), ncol = 1, dimnames = list(rownames(tmat.h[[i]])))
          if(names(tmat.h)[i] == "not_mh") h.flag = FALSE else h.flag = TRUE
          PTmatrix2 = RWR(nadjM = tmat.h[[i]], setSeeds = Seeds, restart = r.param, heterogeneous = h.flag, multiplex = FALSE, net.weight = seed.params$nw)
          rwr_score = as.vector(PTmatrix2)
          # Get ranks of test fold
          ranks = rank(rwr_score)
          tmat.res[[i]][[j]][[l]] = ranks[node.names %in% test.set]
        }
      }
    }
    # Parallel
    if(1){
      cl = makeForkCluster(n.cores, outfile = "")
      registerDoParallel(cl)
      res = foreach(j = 1:length(kegg_cv), .verbose = TRUE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dorng% {
        tmp.res = vector("list", k.folds)
        for(l in seq_len(k.folds)){ # For each fold
          # message(paste0("l",l))
          # Access features of fold l for set j 
          train.set = c(gene_cv[[j]][!gene_fold_ids[[j]] %in% l],
                        comp_cv[[j]][!comp_fold_ids[[j]] %in% l])
          test.set = c(gene_cv[[j]][gene_fold_ids[[j]] %in% l],
                       comp_cv[[j]][comp_fold_ids[[j]] %in% l])
          # Create matrix of seed values
          Seeds = matrix(ifelse(node.names %in% train.set, 1, 0), ncol = 1, dimnames = list(rownames(tmat.h[[i]])))
          if(names(tmat.h)[i] == "not_mh") h.flag = FALSE else h.flag = TRUE
          PTmatrix2 = RWR(nadjM = tmat.h[[i]], setSeeds = Seeds, restart = r.param, heterogeneous = h.flag, multiplex = FALSE, net.weight = seed.params$nw)
          rwr_score = as.vector(PTmatrix2)
          # Get ranks of test fold
          ranks = rank(1 / rwr_score)
          tmp.res[[l]] = ranks[node.names %in% test.set]
        }
        tmp.res
      }
      stopCluster(cl)
      tmat.res[[i]] = res
    }
  }
  saveRDS(tmat.res, file = rank.res.path)
}

#================================#
#=== Monoplex-heterogeneous 2 ===#
#================================#
tmat.res.path = paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/transition_matrix_monoplex_heterogeneous2.rds")
rank.res.path = paste0("/kuhpc/home/s703b506/AMEND2.0/DATA/p", ppi_edge_threshold*10, "m", mmi_edge_threshold*10, "/node.ranking.results_monoplex_heterogeneous2.rds")
if(0){
  message("*** H2 ***")
  # Creating adjacency matrices
  adj.mat.h = igraph::as_adjacency_matrix(g.het2, attr = "weight", sparse = TRUE)
  dimnames(adj.mat.h) = list(paste(rownames(adj.mat.h), V(g.het2)$node_type, sep = "|"),
                             paste(rownames(adj.mat.h), V(g.het2)$node_type, sep = "|"))
  
  # Parameters for RWR
  tmat.params = list(low = list(jp = c(meta = ll, prot = ll)),
                     medium = list(jp = c(meta = mm, prot = mm)),
                     high = list(jp = c(meta = hh, prot = hh)))
  seed.params = list(nw = c(meta = 0.5, prot = 0.5))
  # Construct transition matrices
  if(1){
    tmat.h = vector("list", 4); names(tmat.h) = c("low", "medium", "high", "not_mh")
    for(i in seq_along(tmat.h)){
      message(i)
      if(i == length(tmat.h)){
        tmat.h[[i]] = transition_matrix(adjM = adj.mat.h, norm = "degree", heterogeneous = FALSE, multiplex = FALSE)
      }else{
        tmat.h[[i]] = transition_matrix(adjM = adj.mat.h, norm = "degree", heterogeneous = TRUE, multiplex = FALSE,
                                               jump.prob = tmat.params[[i]]$jp)
      }
    }
    saveRDS(tmat.h, file = tmat.res.path)
  }else tmat.h = readRDS(tmat.res.path)
  
  # Aggregation function
  agg.fun = mean
  # Restart parameter
  r.param = 0.5
  
  #=== k-fold CV ===#
  set.seed(565)
  kegg_cv = kp.h2[unlist(lapply(kp.h2, function(x) length(x$GENE) >= k.folds && length(x$COMPOUND) >= k.folds))]
  gene_cv = lapply(kegg_cv, function(x) unname(x$GENE))
  comp_cv = lapply(kegg_cv, function(x) unname(x$COMPOUND))
  gene_fold_ids = create_folds(kfolds = k.folds, gene_sets = gene_cv)
  comp_fold_ids = create_folds(kfolds = k.folds, gene_sets = comp_cv)
  
  tmat.res = vector("list", length(tmat.h)); names(tmat.res) = names(tmat.h)
  for(i in seq_along(tmat.h)){ # For each transition matrix (low, medium, high, mono-homo)
    message(paste0("i",i))
    node.names = extract_string(rownames(tmat.h[[i]]), "\\|", 1)
    # Not parallel
    if(0){
      tmat.res[[i]] = vector("list", length(kegg_cv))
      for(j in seq_along(kegg_cv)){ # For each pathway
        message(paste0("j",j))
        tmat.res[[i]][[j]] = vector("list", k.folds)
        for(l in seq_len(k.folds)){ # For each fold
          # message(paste0("l",l))
          # Access features of fold l for set j 
          train.set = c(gene_cv[[j]][!gene_fold_ids[[j]] %in% l],
                        comp_cv[[j]][!comp_fold_ids[[j]] %in% l])
          test.set = c(gene_cv[[j]][gene_fold_ids[[j]] %in% l],
                       comp_cv[[j]][comp_fold_ids[[j]] %in% l])
          # Create matrix of seed values
          Seeds = matrix(ifelse(node.names %in% train.set, 1, 0), ncol = 1, dimnames = list(rownames(tmat.h[[i]])))
          if(names(tmat.h)[i] == "not_mh") h.flag = FALSE else h.flag = TRUE
          PTmatrix2 = RWR(nadjM = tmat.h[[i]], setSeeds = Seeds, restart = r.param, heterogeneous = h.flag, multiplex = FALSE, net.weight = seed.params$nw)
          rwr_score = as.vector(PTmatrix2)
          # Get median rank of test fold
          ranks = rank(rwr_score)
          tmat.res[[i]][[j]][[l]] = ranks[node.names %in% test.set]
        }
      }
    }
    # Parallel
    if(1){
      cl = makeForkCluster(n.cores, outfile = "")
      registerDoParallel(cl)
      res = foreach(j = 1:length(kegg_cv), .verbose = TRUE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dorng% {
        tmp.res = vector("list", k.folds)
        for(l in seq_len(k.folds)){ # For each fold
          # message(paste0("l",l))
          # Access features of fold l for set j 
          train.set = c(gene_cv[[j]][!gene_fold_ids[[j]] %in% l],
                        comp_cv[[j]][!comp_fold_ids[[j]] %in% l])
          test.set = c(gene_cv[[j]][gene_fold_ids[[j]] %in% l],
                       comp_cv[[j]][comp_fold_ids[[j]] %in% l])
          # Create matrix of seed values
          Seeds = matrix(ifelse(node.names %in% train.set, 1, 0), ncol = 1, dimnames = list(rownames(tmat.h[[i]])))
          if(names(tmat.h)[i] == "not_mh") h.flag = FALSE else h.flag = TRUE
          PTmatrix2 = RWR(nadjM = tmat.h[[i]], setSeeds = Seeds, restart = r.param, heterogeneous = h.flag, multiplex = FALSE, net.weight = seed.params$nw)
          rwr_score = as.vector(PTmatrix2)
          # Get ranks of test fold
          ranks = rank(1 / rwr_score)
          tmp.res[[l]] = ranks[node.names %in% test.set]
        }
        tmp.res
      }
      stopCluster(cl)
      tmat.res[[i]] = res
    }
  }
  saveRDS(tmat.res, file = rank.res.path)  
}
