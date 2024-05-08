#===================================#
# AMEND Robustness to Degree Bias       
#===================================#
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

library(doParallel)
library(foreach)
# library(doRNG)

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
# AMEND functions
source(paste0(path.amend.code, "run_AMEND.R"))
source(paste0(path.amend.code, "create_integrated_graph.R"))
source(paste0(path.amend.code, "utils.R"))
source(paste0(path.amend.code, "RandomWalk.R"))
source(paste0(path.amend.code, "heinz.R"))

#==========================================#
# Random Degree-Preserving Network Analysis ----
#==========================================#
# RDPN analysis outline:
# For each transition matrix... (5)
#   For each adjustment method... (4)
#     For each GE dataset... (5) 
#       Collect scores on original graph.
#       For each RDPN... (1000) (<-- In Parallel)
#         Collect scores on random graph using permuted seed values.
#       Calculate empirical p-values
# Total number of RWR runs: (1000 + 1) * 5 * 5 * 4 = 100,100
# For figures: display the empirical CDF of p-values. Plot ECDFs of different adj methods on same graph for a given dataset and transition matrix.
#   Perhaps just choose results from one dataset for main figure and include rest as supplemental material.

ll = 0.1
mm = 0.5

n.cores = ifelse(TO_HPC, 8, 3)
r.param = 0.5
data.type = c("pval", "lfc")[1]

path.amend.res = paste0(path.dir, "DB_AMEND_results_", data.type, ".rds")

adj.mat = igraph::as_adjacency_matrix(graph = ppi, attr = "weight", sparse = TRUE)
node.names = rownames(adj.mat)
node.degree = Matrix::colSums(adj.mat)

t.mat.methods = c("degree", paste("modified_degree",ll,sep=";"), paste("modified_degree",mm,sep=";"))
adj.methods = list(control = NULL, SDS = 'SDS', BS = 'BS', IN = 'IN')

# run_AMEND inside parallel for loop
tmat.window = list(1:3)[[1]]
adj.method.window = list(1:4, 1:4, 1:4)

db.res = vector("list", length(t.mat.methods))
names(db.res) = t.mat.methods
for(i in tmat.window){ # seq_along(t.mat.methods)
  if(any(grepl("modified_degree", t.mat.methods[i]))){
    K = as.numeric(extract_string(t.mat.methods[i],";",2))
    tmat.name = "modified_degree"
  }else{
    K = NULL
    tmat.name = t.mat.methods[i]
  } 
  db.res[[i]] = vector("list", length(adj.methods))
  names(db.res[[i]]) = names(adj.methods)
  for(j in adj.method.window[[i]]){ 
    cl = makeForkCluster(n.cores, outfile = "")
    registerDoParallel(cl)
    res = foreach(l = seq_along(omics), .combine = c, .verbose = TRUE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dopar% {
      message(paste0("*** i",i,".j",j,".l",l," ***"))
      v.nm = paste(names(omics)[l],data.type,sep="_")
      seed.values = vertex_attr(ppi, v.nm)
      if(data.type == "pval"){
        seed.values = -log10(seed.values + 1e-6)
      }else if(data.type == "lfc"){
        seed.values = abs(seed.values)
      }
      seed.values[is.na(seed.values)] = 0
      seed.values.asc = sort(seed.values, decreasing = FALSE)
      deg.order.desc = order(node.degree, decreasing = TRUE)
      seed.values[deg.order.desc] = seed.values.asc
      names(seed.values) = node.names
      if(is.null(adj.methods[[j]])) db = NULL else db = list(method = adj.methods[[j]])
      subnet = run_AMEND(adj_matrix = adj.mat, n = 50, data = seed.values, normalize = tmat.name, k = K, degree.bias = db, verbose = TRUE, identifier = paste0("i",i,".j",j,".l",l))
      module.nodes = V(subnet$module)$name
      module.seeds = seed.values[module.nodes]
      module.degree = node.degree[module.nodes]
      module.seeds.ecdf = stats::ecdf(seed.values)(module.seeds)
      module.degree.ecdf = stats::ecdf(node.degree)(module.degree)
      mean(module.seeds.ecdf - module.degree.ecdf)
    }
    stopCluster(cl)
    db.res[[i]][[j]] = res
    names(db.res[[i]][[j]]) = names(omics)
  }
}
saveRDS(db.res, file = paste0(path.results, "DB_AMEND_results.rds"))

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
  
  
  data.type = c("pval", "lfc")[1]
  
  db.res = readRDS(paste0(path.results, "DB_AMEND_results.rds"))

  sum.fun = mean
  val = group = c()
  for(i in seq_along(db.res)){ # tmat
    for(j in seq_along(db.res[[i]])){ # adj method
      tmp = db.res[[i]][[j]]
      val = c(val, sum.fun(tmp)) # take summary 
      tmp1 = ifelse(grepl("modified",names(db.res)[i]), paste0(gsub(";"," (",names(db.res)[i]),")"), names(db.res)[i])
      group = c(group, paste(tmp1,names(db.res[[i]])[j],sep="|"))
    }
  }
  dat = data.frame(val = val,
                   group = group,
                   tmat = factor(extract_string(group, "\\|", 1), levels=c('degree', 'modified_degree (0.1)', 'modified_degree (0.5)')),
                   adjmeth = factor(extract_string(group, "\\|", 2), levels=c('control', 'SDS', 'IN', 'BS')))
  
  col.tmp = brewer.pal(n = length(unique(dat$adjmeth)), name = disc.color.theme)
  
  delta = 30
  plot_title = 'Difference in eCDFs between Seed Valules and Degree of Module Nodes'
  plot_subtitle = paste0('Seed Type: ', ifelse(data.type == 'pval', '-log10(P-Value)', '|LogFC|'))
  y_axis_title = expression(italic(F[s](x) - F[d](x)))
  y.lb = 0
  y.ub = 1
  
  ggplot(dat, aes(x = tmat, y = val, fill = adjmeth)) +
    geom_col(position = "dodge", width=0.7) +
    scale_x_discrete(labels = c("Degree", "Penalized Degree (0.1)", "Penalized Degree (0.5)")) +
    scale_fill_manual(values = col.tmp, labels = c('Control', 'SDS', 'IN', 'BS')) +
    labs(title = plot_title,
         subtitle = plot_subtitle,
         x = "Transition Matrix",
         y = y_axis_title,
         fill = "Adjustment Method") +
    coord_cartesian(ylim=c(y.lb, y.ub)) + 
    theme(plot.title = element_text(size=4+delta,hjust=0,face='bold'),
          plot.subtitle = element_text(size=2+delta,hjust=0,face='bold'),
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






