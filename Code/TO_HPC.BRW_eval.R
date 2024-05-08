#============================#
# Biased Random Walk Analysis
#============================#
# Path names
path.data = "/path/to/data/"
path.amend.code = "/path/to/AMEND/code/"
path.results = "/path/to/results/"

library(igraph)
library(Matrix)
library(doParallel)
library(foreach)
library(doRNG)

# AMEND functions
source(paste0(path.amend.code, "run_AMEND.R"))
source(paste0(path.amend.code, "create_integrated_graph.R"))
source(paste0(path.amend.code, "utils.R"))
source(paste0(path.amend.code, "RandomWalk.R"))
source(paste0(path.amend.code, "heinz.R"))

# Load in DE analysis results (from /Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/PPI Degree Bias/Code/get_Omics.R)
omics = list(als = readRDS(paste0(path.data, "ALS_DE_results.rds")),
             lc = readRDS(paste0(path.data, "LC_DE_results.rds")),
             uc = readRDS(paste0(path.data, "UC_DE_results.rds")),
             cd = readRDS(paste0(path.data, "CD_DE_results.rds")),
             hd = readRDS(paste0(path.data, "HD_DE_results.rds")))
omics = lapply(omics, function(x){
  x$entrez_id = as.numeric(rownames(x))
  x
})

# Load in PPIN (from /Volumes/shared/Biostats/BIO-STAT/Jeffrey Thompson/GRA/Boyd/PPI Degree Bias/Code/get_gene_sets.R)
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

#========================================#
# Active Module Identification with AMEND ----
#========================================#
get_iters = function(id, k){
  a = 1:id
  ca = cut(a, k, labels = FALSE)
  res = numeric(k)
  for(i in seq_len(k)){
    cond = ca %in% unique(ca)[i]
    if(any(cond)) res[i] = min(a[cond]) else res[i] = NA
  } 
  res[!is.na(res)]
}
# get_iters(33, 5)
# get_iters(4,5)

path.orig.mods = "BRW.Eval_original_modules.rds"
path.res.cont = "brw_eval_AMI_cont.rds"
path.res.disc = "brw_eval_AMI_disc.rds"

n.cores = 12

#=== Continuous ===#
# Current
if(0){
  set.seed(12243)
  K = 1:5
  n.cuts = 5
  mod.size = 60
  data.type = c("lfc", "pval")[2]
  
  # Original Modules (no BRW)
  orig.mods = vector("list", length(omics))
  for(i in seq_along(orig.mods)){
    orig.mods[[i]] = run_AMEND(graph = ppi, n = mod.size, data = paste(names(omics)[i], data.type, sep="_"),
                               FUN = "p_value", norm = "degree", brw.attr = NULL, verbose = TRUE, identifier = i)
  }
  saveRDS(orig.mods, file = paste0(path.results, path.orig.mods))
  
  orig.mods = readRDS(paste0(path.results, path.orig.mods))
  
  brw.res = vector("list", n.cuts)
  for(i in seq_len(n.cuts)){ # 5
    message(paste0("i",i))
    brw.res[[i]] = vector("list", length(K))
    for(j in seq_along(K)){ # 5
      message(paste0("j",j))
      brw.res[[i]][[j]] = vector("list", 3); names(brw.res[[i]][[j]]) = c("Iter_rm", "Dist", "N")
      brw.res[[i]][[j]] = lapply(brw.res[[i]][[j]], function(x) numeric(length(omics)))
      # Not parallel
      if(0){
        for(l in seq_along(omics)){ # 5
          message(paste0("l",l))
          best.id = which.max(orig.mods[[l]]$stats$`Network score`)
          iter = get_iters(id = best.id, k = n.cuts)
          iter = iter[min(length(iter), i)]
          noi = orig.mods[[l]]$subnetworks[[iter]]
          # noi = noi[!noi %in% V(orig.mods[[l]]$module)$name]
          noi = noi[!noi %in% orig.mods[[l]]$subnetworks[[iter+1]]]
          ids = which(V(ppi)$name %in% noi)
          brw.values = numeric(vcount(ppi))
          brw.values[ids] = rexp(length(ids), 1/K[j]) + 1
          brw.values[-ids] = 1
          names(brw.values) = V(ppi)$name
          brw.mod = run_AMEND(graph = ppi, n = mod.size, data = paste(names(omics)[l], data.type, sep="_"),
                              FUN = "p_value", norm = "degree", brw.attr = brw.values, verbose = TRUE, identifier = l)
          nn = V(brw.mod$module)$name
          best.id.brw = which.max(brw.mod$stats$`Network score`)
          # Average increase in iteration removal from orig to brw modules
          # Scale/Normalize diff.iter.rm by the iteration number of final modules
          diff.iter.rm = numeric(length(nn))
          for(ii in seq_along(nn)){
            tmp.orig = max(which(unlist(lapply(orig.mods[[l]]$subnetworks, function(y) nn[ii] %in% y )))) 
            tmp.brw = max(which(unlist(lapply(brw.mod$subnetworks, function(y) nn[ii] %in% y )))) 
            diff.iter.rm[ii] = tmp.brw / best.id.brw - tmp.orig / best.id
          }
          diff.iter.rm = mean(diff.iter.rm)
          # Average decrease in average distance between NOIs and module nodes from orig to brw modules
          d.mat.orig = igraph::distances(graph = ppi, v = noi, to = V(orig.mods[[l]]$module)$name)
          d.mat.brw = igraph::distances(graph = ppi, v = noi, to = V(brw.mod$module)$name)
          diff.dist = mean(apply(d.mat.orig, 1, mean) - apply(d.mat.brw, 1, mean))
          brw.res[[i]][[j]][["Iter_rm"]][l] = diff.iter.rm
          brw.res[[i]][[j]][["Dist"]][l] = diff.dist
          brw.res[[i]][[j]][["N"]][l] = length(noi)
        }
      }
      # Parallel
      if(1){
        cl = makeForkCluster(n.cores, outfile = "")
        registerDoParallel(cl)
        res = foreach(l = 1:length(omics), .verbose = TRUE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dorng% {
          message(paste0("l",l))
          best.id = which.max(orig.mods[[l]]$stats$`Network score`)
          iter = get_iters(id = best.id, k = n.cuts)
          iter = iter[min(length(iter), i)]
          noi = orig.mods[[l]]$subnetworks[[iter]]
          # noi = noi[!noi %in% V(orig.mods[[l]]$module)$name]
          noi = noi[!noi %in% orig.mods[[l]]$subnetworks[[iter+1]]]
          ids = which(V(ppi)$name %in% noi)
          brw.values = numeric(vcount(ppi))
          brw.values[ids] = rexp(length(ids), 1/K[j]) + 1
          brw.values[-ids] = 1
          names(brw.values) = V(ppi)$name
          brw.mod = run_AMEND(graph = ppi, n = mod.size, data = paste(names(omics)[l], data.type, sep="_"),
                              FUN = "p_value", norm = "degree", brw.attr = brw.values, verbose = TRUE, identifier = l)
          nn = V(brw.mod$module)$name
          best.id.brw = which.max(brw.mod$stats$`Network score`)
          # Average increase in iteration removal from orig to brw modules
          # Scale/Normalize diff.iter.rm by the iteration number of final modules
          diff.iter.rm = numeric(length(nn))
          for(ii in seq_along(nn)){
            tmp.orig = max(which(unlist(lapply(orig.mods[[l]]$subnetworks, function(y) nn[ii] %in% y )))) 
            tmp.brw = max(which(unlist(lapply(brw.mod$subnetworks, function(y) nn[ii] %in% y )))) 
            diff.iter.rm[ii] = tmp.brw / best.id.brw - tmp.orig / best.id
          }
          diff.iter.rm = mean(diff.iter.rm)
          # Average decrease in average distance between NOIs and module nodes from orig to brw modules
          d.mat.orig = igraph::distances(graph = ppi, v = noi, to = V(orig.mods[[l]]$module)$name)
          d.mat.brw = igraph::distances(graph = ppi, v = noi, to = V(brw.mod$module)$name)
          diff.dist = mean(apply(d.mat.orig, 1, mean) - apply(d.mat.brw, 1, mean))
          tmp = c(diff.iter.rm, diff.dist, length(noi))
          names(tmp) = c("Iter_rm", "Dist", "N")
          tmp
        }
        stopCluster(cl)
        # change result format
        for(l in seq_along(omics)) 
          for(ii in seq_along(brw.res[[i]][[j]])) 
            brw.res[[i]][[j]][[ii]][l] = res[[l]][ii]
      }
    }
  }
  saveRDS(brw.res, file = paste0(path.results, path.res.cont))
}
# bar plot
if(0){
  if(0){
    display.brewer.all()
    display.brewer.pal(n = 5, name = div.color.theme)
    brewer.pal(n = 5, name = div.color.theme)
  }
  
  lags = 1:5 # from beginning to end 
  inv.rates = 1:5
  .iter = .dist = group = size = c()
  sum.fun = mean
  for(i in seq_along(brw.res.cont)){
    for(j in seq_along(brw.res.cont[[i]])){
      tmp1 = unlist(brw.res.cont[[i]][[j]][[1]])
      .iter = c(.iter, sum.fun(tmp1))
      tmp2 = unlist(brw.res.cont[[i]][[j]][[2]])
      .dist = c(.dist, sum.fun(tmp2))
      group = c(group, paste(lags[i], inv.rates[j], sep="_"))
      size = c(size, median(brw.res.cont[[i]][[j]][["N"]]))
    }
  }
  dat = data.frame(iter = .iter,
                   dist = -.dist,
                   lag = extract_string(group,"_",1),
                   inv.rate = extract_string(group,"_",2),
                   N = size,
                   row.names=NULL)
  
  # Plot Settings
  col.tmp = brewer.pal(n = length(unique(dat$inv.rate)), name = disc.color.theme)
  plot_title = c("AMEND: Change in Iteration Removal Number","")[1]
  plot_subtitle = c("from RWR to B-RWR","")[1]
  y_label = "\u0394-Iteration-of-Removal (scaled)"
  delta = 20
  
  p.iter = ggplot(dat, aes(x = lag, y = iter, fill = inv.rate)) +
    geom_col(position = "dodge", width=0.51) +
    theme_light() +
    scale_fill_manual(values = col.tmp, labels = c("1/1", "1/2", "1/3", "1/4", "1/5")) +
    theme(plot.title = element_text(size=12+delta,face='bold'),
          plot.subtitle = element_text(size=11+delta,face='bold'),
          legend.text = element_text(size=10+delta,face='bold.italic'),
          legend.title = element_text(size=10+delta,face='bold'),
          axis.text.y = element_text(size=10+delta,face='bold.italic'),
          axis.title.y = element_text(size=12+delta,face='bold',margin=margin(r=10)),
          axis.title.x = element_text(size=12+delta,face='bold',margin=margin(t=10)),
          axis.ticks.x = element_line(lineend=c('round','butt','square')[3],color='black',linewidth = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last",type='closed')),
          axis.text.x = element_text(size=12+delta,face='bold.italic')) +
    labs(title = plot_title,
         subtitle = plot_subtitle,
         x = "AMEND Iterations",
         y = y_label,
         fill = "Rate") +
    scale_x_discrete(labels = c("Start", "", "", "", "End"))
  p.iter
  
  ###
  plot_title = c("Mean \u0394-Distance for Node Sets", 
                 "Mean Decrease in Distance from Node Set to Module Nodes",
                 "AMEND: Change in Distance between NOI and Module Nodes")[3]
  plot_subtitle = c("Original Module vs. BRW Module",
                    "from RWR to B-RWR")[2]
  y_label = "\u0394-Distance"
  
  p.dist = ggplot(dat, aes(x = lag, y = dist, fill = inv.rate)) +
    geom_col(position = "dodge", width=0.51) +
    theme_light() +
    scale_fill_manual(values = col.tmp, labels = c("1/1", "1/2", "1/3", "1/4", "1/5")) +
    theme(plot.title = element_text(size=12+delta,face='bold'),
          plot.subtitle = element_text(size=11+delta,face='bold'),
          legend.text = element_text(size=10+delta,face='bold.italic'),
          legend.title = element_text(size=10+delta,face='bold'),
          axis.text.y = element_text(size=10+delta,face='bold.italic'),
          axis.title.y = element_text(size=12+delta,face='bold',margin=margin(l=10)),
          axis.title.x = element_text(size=12+delta,face='bold',margin=margin(t=10)),
          axis.ticks.x = element_line(lineend=c('round','butt','square')[3],color='black',linewidth = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last",type='closed')),
          axis.text.x = element_text(size=12+delta,face='bold.italic')) +
    labs(title = plot_title,
         subtitle = plot_subtitle,
         x = "AMEND Iterations",
         y = y_label,
         fill = "Rate") +
    scale_x_discrete(labels = c("Start", "", "", "", "End")) +
    geom_hline(yintercept=0,color='lightgrey',linewidth=0.8)
  p.dist
  
  # Save dimensions: 
  # width=900, height=700
  # size factor of 2 --. width=1800, height=1400, with delta=20
}

#=== Discrete ===#
# Current
if(0){
  n.cuts = 5
  mod.size = 60
  data.type = c("lfc", "pval")[2]
  
  # Original Modules (no BRW)
  orig.mods = readRDS(path.orig.mods)
  
  brw.res = vector("list", n.cuts)
  for(i in seq_len(n.cuts)){ # 5
    message(paste0("i",i))
    brw.res[[i]] = vector("list", 3); names(brw.res[[i]]) = c("Iter_rm", "Dist", "N")
    brw.res[[i]] = lapply(brw.res[[i]], function(x) numeric(length(omics)))
    # Not parallel
    if(0){
      for(l in seq_along(omics)){ # 5
        message(paste0("l",l))
        best.id = which.max(orig.mods[[l]]$stats$`Network score`)
        iter = get_iters(id = best.id, k = n.cuts)
        iter = iter[min(length(iter), i)]
        noi = orig.mods[[l]]$subnetworks[[iter]]
        noi = noi[!noi %in% orig.mods[[l]]$subnetworks[[iter+1]]]
        brw.mod = run_AMEND(graph = ppi, n = mod.size, data = paste(names(omics)[l], data.type, sep="_"),
                            FUN = "p_value", norm = "degree", brw.attr = noi, verbose = TRUE, identifier = l)
        nn = V(brw.mod$module)$name
        best.id.brw = which.max(brw.mod$stats$`Network score`)
        # Average increase in iteration removal from orig to brw modules
        diff.iter.rm = numeric(length(nn))
        for(ii in seq_along(nn)){
          tmp.orig = max(which(unlist(lapply(orig.mods[[l]]$subnetworks, function(y) nn[ii] %in% y ))))
          tmp.brw = max(which(unlist(lapply(brw.mod$subnetworks, function(y) nn[ii] %in% y ))))
          diff.iter.rm[ii] = tmp.brw / best.id.brw - tmp.orig / best.id
        }
        diff.iter.rm = mean(diff.iter.rm)
        # Average decrease in average distance between NOIs and module nodes from orig to brw modules
        d.mat.orig = igraph::distances(graph = ppi, v = noi, to = V(orig.mods[[l]]$module)$name)
        d.mat.brw = igraph::distances(graph = ppi, v = noi, to = V(brw.mod$module)$name)
        diff.dist = mean(apply(d.mat.orig, 1, mean) - apply(d.mat.brw, 1, mean))
        brw.res[[i]][["Iter_rm"]][l] = diff.iter.rm
        brw.res[[i]][["Dist"]][l] = diff.dist
        brw.res[[i]][["N"]][l] = length(noi)
      }
    }
    # Parallel
    if(1){
      cl = makeForkCluster(n.cores, outfile = "")
      registerDoParallel(cl)
      res = foreach(l = 1:length(omics), .verbose = TRUE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dorng% {
        message(paste0("l",l))
        best.id = which.max(orig.mods[[l]]$stats$`Network score`)
        iter = get_iters(id = best.id, k = n.cuts)
        iter = iter[min(length(iter), i)]
        noi = orig.mods[[l]]$subnetworks[[iter]]
        noi = noi[!noi %in% orig.mods[[l]]$subnetworks[[iter+1]]]
        brw.mod = run_AMEND(graph = ppi, n = mod.size, data = paste(names(omics)[l], data.type, sep="_"),
                            FUN = "p_value", norm = "degree", brw.attr = noi, verbose = TRUE, identifier = l)
        nn = V(brw.mod$module)$name
        best.id.brw = which.max(brw.mod$stats$`Network score`)
        # Average increase in iteration removal from orig to brw modules
        diff.iter.rm = numeric(length(nn))
        for(ii in seq_along(nn)){
          tmp.orig = max(which(unlist(lapply(orig.mods[[l]]$subnetworks, function(y) nn[ii] %in% y ))))
          tmp.brw = max(which(unlist(lapply(brw.mod$subnetworks, function(y) nn[ii] %in% y ))))
          diff.iter.rm[ii] = tmp.brw / best.id.brw - tmp.orig / best.id
        }
        diff.iter.rm = mean(diff.iter.rm)
        # Average decrease in average distance between NOIs and module nodes from orig to brw modules
        d.mat.orig = igraph::distances(graph = ppi, v = noi, to = V(orig.mods[[l]]$module)$name)
        d.mat.brw = igraph::distances(graph = ppi, v = noi, to = V(brw.mod$module)$name)
        diff.dist = mean(apply(d.mat.orig, 1, mean) - apply(d.mat.brw, 1, mean))
        tmp = c(diff.iter.rm, diff.dist, length(noi))
        names(tmp) = c("Iter_rm", "Dist", "N")
        tmp
      }
      stopCluster(cl)
      # change result format
      for(l in seq_along(omics)) 
        for(j in seq_along(brw.res[[i]])) 
          brw.res[[i]][[j]][l] = res[[l]][j] 
    }
  }
  saveRDS(brw.res, file = paste0(path.results, path.res.disc))
}
# bar plot
if(0){
  lags = 1:5 # from beginning to end 
  .iter = .dist = group = size = c()
  sum.fun = mean
  for(i in seq_along(brw.res.disc)){
    tmp1 = unlist(brw.res.disc[[i]][[1]])
    .iter = c(.iter, sum.fun(tmp1))
    tmp2 = unlist(brw.res.disc[[i]][[2]])
    .dist = c(.dist, sum.fun(tmp2))
    group = c(group, lags[i])
    size = c(size, median(brw.res.disc[[i]][["N"]]))
  }
  dat = data.frame(iter = .iter,
                   dist = -.dist,
                   lag = group,
                   N = size,
                   row.names=NULL)
  
  # Distance
  plot_title = c("AMEND: Change in Distance between NOI and Module Nodes")[1]
  plot_subtitle = c("from RWR to B-RWR")[1]
  y_label = "\u0394-Distance"
  delta = 20
  
  p.dist = ggplot(dat, aes(x = factor(lag), y = dist)) +
    geom_col(position = "dodge", width=0.51, fill = 'blue3') +
    geom_hline(yintercept=0,color='darkgrey',linewidth=0.8) +
    scale_x_discrete(labels = c("Start", "", "", "", "End")) +
    theme_light() +
    theme(plot.title=element_text(size=12+delta,face='bold'),
          plot.subtitle=element_text(size=11+delta,face='bold'),
          axis.title.y=element_text(size=12+delta,face='bold',margin=margin(r=10)),
          axis.text.y = element_text(size=10+delta,face='bold.italic'),
          axis.title.x = element_text(size=12+delta,face='bold',margin=margin(t=10)),
          axis.ticks.x = element_line(lineend=c('round','butt','square')[3],color='black',linewidth = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last",type='closed')),
          axis.text.x = element_text(size=12+delta,face='bold.italic')) +
    labs(title = plot_title,
         subtitle = plot_subtitle,
         x = "AMEND Iterations",
         y = y_label) 
  p.dist
  
  # Iteration removal
  plot_title = c("AMEND: Change in Iteration Removal Number")[1]
  plot_subtitle = c("from RWR to B-RWR")[1]
  y_label = "\u0394-Iteration-of-Removal (scaled)"
  
  p.iter = ggplot(dat, aes(x = factor(lag), y = iter)) +
    geom_col(position = "dodge", width=0.51, fill = 'blue3') +
    geom_hline(yintercept=0,color='darkgrey',linewidth=0.8) +
    scale_x_discrete(labels = c("Start", "", "", "", "End")) +
    theme_light() +
    theme(plot.title=element_text(size=12+delta,face='bold'),
          plot.subtitle=element_text(size=11+delta,face='bold'),
          axis.title.y=element_text(size=12+delta,face='bold',margin=margin(r=10)),
          axis.text.y = element_text(size=10+delta,face='bold.italic'),
          axis.title.x = element_blank(),
          axis.ticks.x = element_line(lineend=c('round','butt','square')[3],color='black',linewidth = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last",type='closed')),
          axis.text.x = element_text(size=12+delta,face='bold.italic')) +
    labs(title = plot_title,
         subtitle = plot_subtitle,
         x = "AMEND Iterations",
         y = y_label) 
  p.iter
  
  plots = plot_grid(plotlist = list(p.iter, p.dist), align = "v", nrow = 2)
  plots
  
  # Save dimensions: 
  # width=900, height=700
  # size factor of 2 --. width=1800, height=1400, with delta=20
}

#===============#
# RWR ----
#===============#
if(0){
  set.seed(125035)
  from_clusters = TRUE
  data.type = c("pval", "lfc")[1] 
  r.param = 0.5 
  
  #=== Louvain Clustering ===#
  clusts = igraph::cluster_louvain(ppi)
  clusts$membership[1:5]
  c.ids = unique(clusts$membership)
  good.ids = c.ids[c.ids %in% names(table(clusts$membership))[table(clusts$membership) >= n]]
  
  ppi.adj = igraph::as_adjacency_matrix(ppi, attr = "weight", sparse = TRUE)
  brw.res = vector("list", length(omics))
  for(i in seq_along(omics)){ # 5
    message(paste0("i",i))
    # RWR scores without BRW
    tmat.orig = transition_matrix(adjM = ppi.adj, norm = "degree")
    seeds = matrix(igraph::vertex_attr(ppi, paste(names(omics)[i], data.type, sep="_")),ncol=1,dimnames=list(V(ppi)$name))
    if(data.type == "lfc") seeds = abs(seeds) else seeds = -log10(seeds+0.01)
    seeds[is.na(seeds) | seeds < 0] = 0
    rwr.orig = RWR(nadjM = tmat.orig, setSeeds = seeds, restart = r.param)[,1]
    rwr.orig.quant = stats::ecdf(rwr.orig)(rwr.orig)
    rwr.orig.ranks = rank(1 / rwr.orig)
    brw.res[[i]] = vector("list", length(K))
    for(j in seq_along(K)){ # 5
      message(paste0("j",j))
      res = numeric(B)
      for(l in seq_len(B)){ # 100
        randnums = rexp(n, 1/K[j]) + 1
        if(!from_clusters){
          ids = sample(1:N, n)
        }else{
          c.id.tmp = sample(good.ids, 1)
          ids = sample(which(clusts$membership == c.id.tmp), n)
        }
        brw.values = numeric(N)
        brw.values[ids] = randnums
        brw.values[-ids] = 1
        tmat = transition_matrix(adjM = ppi.adj, norm = "degree", brw.attr = brw.values)
        # if(!all(rownames(tmat) == rownames(ppi.adj))) stop("problem")
        rwr.brw = RWR(nadjM = tmat, setSeeds = seeds, restart = r.param)[,1]
        rwr.brw.quant = stats::ecdf(rwr.brw)(rwr.brw)
        # rwr.brw.ranks = rank(1 / rwr.brw) 
        # res[l] = median(rwr.brw[ids] - rwr.orig[ids])
        # res[l] = median(rwr.orig.ranks[ids] - rwr.brw.ranks[ids]) # not scaled
        # res[l] = median(rwr.orig.ranks[ids] - rwr.brw.ranks[ids]) / N # scaled by size of graph
        res[l] = median(rwr.brw.quant[ids] - rwr.orig.quant[ids])
      }
      brw.res[[i]][[j]] = res
    }
  }
  names(brw.res) = names(omics)
  brw.res = lapply(brw.res, function(x) {
    names(x) = K
    x
  })
  saveRDS(brw.res, file = paste0(path.results, "brw.rwr.res.rds"))
}

# Line plot
if(0){
  extract_string = function(x,k,pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))
  agg.rank = v.names = min.rank = max.rank = p25.rank = p75.rank = NULL
  sum.fun = median
  for(i in seq_along(brw.res)){
    for(j in seq_along(brw.res[[i]])){
      agg.rank = c(agg.rank, sum.fun(brw.res[[i]][[j]]))
      min.rank = c(min.rank, min(brw.res[[i]][[j]]))
      max.rank = c(max.rank, max(brw.res[[i]][[j]]))
      p25.rank = c(p25.rank, quantile(brw.res[[i]][[j]], 0.25))
      p75.rank = c(p75.rank, quantile(brw.res[[i]][[j]], 0.75))
      nm = paste(names(brw.res)[i], names(brw.res[[i]])[j], sep="_")
      v.names = c(v.names, nm)
    }
  }
  dat = data.frame(agg.rank = agg.rank,
                   data = extract_string(v.names,"_",1),
                   rate = as.numeric(extract_string(v.names,"_",2)),
                   row.names = NULL)
  
  plot_title = "Change in eCDF"
  plot_subtitle = "from RWR to B-RWR"
  y_label = "\u0394-eCDF"
  bar.colors = brewer.pal(n = 5, name = disc.color.theme)
  delta = 20
  
  brw.plot1 = ggplot(dat, aes(rate, agg.rank, color = data)) +
    geom_line(linewidth = 1) +
    geom_point(shape = "diamond", size = 5) +
    labs(title = plot_title, subtitle = plot_subtitle, x = "Exponential Rate", y = y_label, color = "Dataset") +
    scale_color_manual(values = bar.colors, labels = toupper(sort(unique(dat$data)))) +
    scale_x_continuous(labels = c("1/1", "1/2", "1/3", "1/4", "1/5")) +
    theme(plot.title = element_text(hjust=0,size=12+delta,face="bold"),
          plot.subtitle = element_text(hjust=0,size=11+delta,face='bold'),
          axis.title.x.bottom = element_text(size = 12+delta, face = "bold"),
          axis.title.y.left = element_text(size = 12+delta, face = "bold"),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major.y = element_line(colour = "black", linetype = 2, linewidth = 0.2),
          axis.text = element_text(size = 12+delta, face = "bold.italic"),
          legend.title = element_text(size = 12+delta, face = "bold"), 
          legend.text = element_text(size = 8+delta, face = "bold"))
  
  brw.plot1
  # Save dimensions: 
  # width=900, height=700
  # size factor of 2 --. width=1800, height=1400, with delta=20
}
