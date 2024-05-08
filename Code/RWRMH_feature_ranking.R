#================#
# RWR-MH vs. RWR
# Feature Ranking
#================#

# To-Do:
#   - Construct a MH graph... since I'm introducing RWR-MH, I should analyze a MH graph, not just a heterogeneous graph.
#   - Use transition_matrix() with M=TRUE,H=TRUE, M=FALSE,H=TRUE, or M=FALSE,H=FALSE
#   - Create an 'agg_name' vertex attribute, then use code in 'restart_grid_search'
#   - Write code to get average rank of nodes in hold-out set for each pathway, for both MH and H graphs

# Clear environment
rm(list = ls())

# Path names
path.data = "/path/to/data/"
path.amend.code = "/path/to/AMEND/code/"
path.results = "/path/to/results/"

# Attach libraries
library(KEGGREST) # R to KEGG API
library(data.table) # Fast import of files
library(dplyr) # Easy data manipulation
library(igraph) # graph creation and manipulation
library(biomaRt) # Annotation mapping
library(graphite) # Pathway network
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(RColorBrewer)

# AMEND functions
source(paste0(path.amend.conde, "run_AMEND.R"))
source(paste0(path.amend.conde, "RandomWalk.R"))
source(paste0(path.amend.conde, "create_integrated_graph.R"))
source(paste0(path.amend.conde, "utils.R"))

if(0){
  display.brewer.all()
  display.brewer.pal(n = 5, name = "Set1")
  brewer.pal(n = 5, name = "Set1")
}

disc.color.theme = "Set1"
div.color.theme = "PRGn"
grad.color.theme ="Greens"

# Functions
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

#============================#
# Get KEGG Metabolic Pathways ----
#============================#
pathnames = keggList("pathway", "mmu")
res = vector("list", length(pathnames))
kegg_pathways = list()
for(i in seq_along(res)){
  if(i %% 50 == 0) message(i)
  tmp = keggGet(names(pathnames)[i])
  res[[i]] = names(tmp[[1]])
  if(all(c("COMPOUND", "GENE") %in% res[[i]])) kegg_pathways[[length(kegg_pathways)+1]] = tmp[[1]]
}
table(unlist(res))

kegg_pathways[[1]]$GENE[1:10]
kegg_pathways[[1]]$COMPOUND[1:10]

# For genes, set Entrez ID as name, set gene symbol as value
for(i in seq_along(kegg_pathways)){
  genes = kegg_pathways[[i]]$GENE
  n = length(genes)
  if(n %% 2 != 0) stop("uneven length of gene set.")
  entrez_id = genes[(1:(n/2))*2-1]
  gene_symbol = extract_string(genes[(1:(n/2))*2], ";", 1)
  kegg_pathways[[i]]$GENE = gene_symbol
  names(kegg_pathways[[i]]$GENE) = entrez_id
}
# saveRDS(kegg_pathways, file = "/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Data/kegg_pathways.rds")
# Genes as Entrez IDs and Gene Symbols. Metabolites as KEGG Compound IDs and Compound Names
# kegg_pathways = readRDS(file = "/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Data/kegg_pathways.rds")

#=================================================#
# Construct Protein-Metabolite Interaction Network ----
#=================================================#
# Large network files will be stored in S-drive
ppi_edge_threshold = 0.9
mmi_edge_threshold = 0.9

#=== Physical PPI Network ===#
ppi = data.table::fread(file = paste0(path.data, "10090.protein.physical.links.v12.0.txt.gz"), header = T) %>%
  dplyr::mutate(combined_score = combined_score / 1000) %>%
  dplyr::filter(combined_score >= ppi_edge_threshold) %>%
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
# Save
# saveRDS(ppi, file = paste0(path.results, "ppi_edgelist_", ppi_edge_threshold, ".rds"))

#=== MMI Network of merged isomers ===#
# CIDs = stereo-specific isomer (i.e., spatially different isomers)
# CIDm = 'flat' compounds (i.e., merged stereo-isomers)
mmi = data.table::fread(file = paste0(path.data, "chemical_chemical.links.detailed.v5.0.tsv.gz"), header = T, sep = "\t") %>%
  dplyr::mutate(combined_score = combined_score / 1000) %>%
  dplyr::filter(combined_score >= mmi_edge_threshold) %>%
  dplyr::filter(substr(chemical1,1,4) == "CIDm" & substr(chemical2,1,4) == "CIDm") %>% # Keep only merged isomers
  dplyr::mutate(chemical1 = gsub(pattern = "^0+", replacement = "", x = substr(chemical1, 5, 555)), # remove CIDm/CIDs and leading zeros
                chemical2 = gsub(pattern = "^0+", replacement = "", x = substr(chemical2, 5, 555))) %>%
  dplyr::filter(chemical1 != chemical2) %>%
  dplyr::select(chemical1, chemical2, combined_score) %>%
  dplyr::distinct() %>%
  as.matrix()

# Remove duplicates and Aggregate edge scores
tmp = character(nrow(mmi))
for(i in seq_along(tmp)){
  tmp[i] = paste(sort(c(mmi[i,1], mmi[i,2])), collapse = "|")
}
scores = as.numeric(mmi[,3])
dat = stats::aggregate(scores, by = list(tmp), FUN = min)
mmi = matrix(c(extract_string(dat[,1], "\\|", 1), extract_string(dat[,1], "\\|", 2), dat[,2]), ncol = 3)
# Save
# saveRDS(mmi, file = paste0(path.results, "mmi_edgelist_", mmi_edge_threshold, ".rds"))

#=== Protein-Metabolite Interactions ===#
bp = data.table::fread(file = paste0(path.data, "10090.protein_chemical.links.detailed.v5.0.tsv.gz"), header = T, sep = "\t") %>%
  dplyr::mutate(combined_score = combined_score / 1000) %>%
  dplyr::filter(combined_score >= ppi_edge_threshold) %>%
  dplyr::filter(substr(chemical,1,4) == "CIDm") %>% # Keep only merged isomers
  dplyr::mutate(chemical = gsub(pattern = "^0+", replacement = "", x = substr(chemical, 5, 555)),
                protein = extract_string(protein, "\\.", 2)) %>%
  dplyr::select(chemical, protein, combined_score) %>%
  dplyr::distinct() %>%
  as.matrix()
colnames(bp) = NULL
# Save
# saveRDS(bp, file = paste0(path.results, "bipartite_edgelist_", ppi_edge_threshold, ".rds"))

#=== Pathway Network from Graphite ===#
graphite::pathwayDatabases()
mmR = graphite::pathways("mmusculus", "reactome")
for(i in seq_along(mmR)){ 
  if(i %% 100 == 0) message(i)
  g.tmp = graphite::pathwayGraph(pathway = mmR[[i]], which = "proteins")
  g.tmp = igraph::graph_from_graphnel(graphNEL = g.tmp)
  g.tmp = igraph::as.undirected(graph = g.tmp, mode = "collapse")
  if(i == 1) g.path = g.tmp else g.path = igraph::union(g.path, g.tmp)
}
edge.weight.attrs = igraph::edge_attr_names(g.path)[grepl("weight", igraph::edge_attr_names(g.path))]
for(i in seq_along(edge.weight.attrs)){
  g.path = igraph::delete_edge_attr(graph = g.path, name = edge.weight.attrs[i])
}
E(g.path)$weight = ppi_edge_threshold # Give all edges the minimum edge weight
V(g.path)$name = gsub("UNIPROT:", "", V(g.path)$name)
uniq.uniprot = unique(V(g.path)$name)
# Map UNIPROT IDs to Gene Symbols
if(1){
  mart <- useMart(biomart = "ensembl")
  mart_data <- useDataset("mmusculus_gene_ensembl", mart = mart)
  # Check the names of available filters & attributes
  if(0){
    View(listFilters(mart_data))
    View(listAttributes(mart_data))  
  }
  # Attributes are what you want returned from your mapping (i.e., what you are mapping to)
  # Filters are what you are trying to map 
  # Example: If you want to map Symbols to Entrez_IDs, then filters=Symbols, attributes=c(Symbols, Entrez_IDs). attributes has both because we want to return both types.
  gene.map <- getBM(attributes = c("mgi_symbol", "uniprot_gn_id"),
                    filters = "uniprot_gn_id",
                    values = uniq.uniprot,
                    mart = mart_data)
  gene.map = gene.map[gene.map[,1] != "" & gene.map[,2] != "", ] # Remove features that didn't map
  head(gene.map); nrow(gene.map)
  # Check relationship of mapping (e.g., many-to-one, many-to-many, etc.)
  table(table(gene.map[,1]))
  table(table(gene.map[,2]))
  # Remove UNIPROT IDs from gene.map that have multiple matches
  multi.uniprot = names(table(gene.map[,2]))[table(gene.map[,2]) > 1]
  gene.map = gene.map[!gene.map[,2] %in% multi.uniprot,]
  # Remove Nodes in g.path that didn't map
  g.path = igraph::induced_subgraph(graph = g.path, vids = which(V(g.path)$name %in% gene.map[,2]))
  # Mapping symbols to graph
  V(g.path)$uniprot = V(g.path)$name
  ids = match(V(g.path)$uniprot, gene.map[,2])
  igraph::vertex_attr(g.path, "name", which(!is.na(ids))) = gene.map[ids[!is.na(ids)],1]
}
# saveRDS(g.path, file = "pathway_graph.rds")
# g.path = readRDS("pathway_graph.rds")
E(g.path)$weight = ppi_edge_threshold # Give all edges the minimum edge weight

#=== Heterogeneous Network: Physical PPI ----
uniq.prot.phys = unique(c(ppi[,1], ppi[,2]))
uniq.meta = unique(c(mmi[,1], mmi[,2]))
bp1 = bp[bp[,1] %in% uniq.meta & bp[,2] %in% uniq.prot.phys,]

het_el = rbind(ppi, mmi, bp1)

g.het1 = igraph::graph_from_edgelist(el = het_el[,1:2], directed = FALSE)
E(g.het1)$weight = as.numeric(het_el[,3])
V(g.het1)$node_type = ifelse(grepl("ENSMUSP", V(g.het1)$name), "prot", "meta")

# Unique compound names and KEGG IDs in KEGG pathways
uniq.kegg.comp = unique(unlist(lapply(kegg_pathways, function(x) unname(x$COMPOUND))))
uniq.kegg.id = unique(unlist(lapply(kegg_pathways, function(x) names(x$COMPOUND))))

# Unique Entrez IDs and Gene symbols in KEGG pathways
uniq.kegg.entrez = unique(unlist(lapply(kegg_pathways, function(x) names(x$GENE))))
uniq.kegg.symbol = unique(unlist(lapply(kegg_pathways, function(x) unname(x$GENE))))

#=== Mapping with biomaRt ===#
## Proteins/Genes
# ENSMUSP IDs to Gene Symbols
if(1){
  mart <- useMart(biomart = "ensembl")
  mart_data <- useDataset("mmusculus_gene_ensembl", mart = mart)
  # Check the names of available filters & attributes
  if(0){
    View(listFilters(mart_data))
    View(listAttributes(mart_data))  
  }
  # Attributes are what you want returned from your mapping (i.e., what you are mapping to)
  # Filters are what you are trying to map 
  # Example: If you want to map Symbols to Entrez_IDs, then filters=Symbols, attributes=c(Symbols, Entrez_IDs). attributes has both because we want to return both types.
  prot.map <- getBM(attributes = c("ensembl_peptide_id", "mgi_symbol"),
                    filters = "ensembl_peptide_id",
                    values = uniq.prot.phys,
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
  g.het1 = igraph::induced_subgraph(g.het1, which(V(g.het1)$name %in% prot.map[,1] | V(g.het1)$node_type == "meta"))
  # Map Symbols to nodes in network
  id = match(V(g.het1)$name, prot.map[,1])
  igraph::vertex_attr(g.het1, "mgi_symbol", which(!is.na(id))) = prot.map[id[!is.na(id)],2]
}

## Compounds
# Compound Name to PubChemID 
comp.map = data.table::fread(paste0(path.data, "MetaboAnalyst_Compound_Name.txt"), sep = ",")

igraph::vertex_attr(g.het1, "comp_name", which(V(g.het1)$node_type == "meta")) = V(g.het1)$name[V(g.het1)$node_type == "meta"]
igraph::vertex_attr(g.het1, "ensp", which(V(g.het1)$node_type == "prot")) = V(g.het1)$name[V(g.het1)$node_type == "prot"]
igraph::vertex_attr(g.het1, "name", which(V(g.het1)$node_type == "prot")) = V(g.het1)$mgi_symbol[V(g.het1)$node_type == "prot"]

# Take largest connected component
g.het1 = largest_connected_component(g.het1)

V(g.het1)$label = paste(V(g.het1)$name, V(g.het1)$node_type, sep = "|")

#=== Heterogeneous Network: Functional PPI ----
bp2 = bp[bp[,1] %in% uniq.meta,]
uniq.prot.bp = unique(bp2[,2])

#=== Mapping ===#
# ENSMUSP IDs to Gene Symbols
mart <- useMart(biomart = "ensembl")
mart_data <- useDataset("mmusculus_gene_ensembl", mart = mart)
# Check the names of available filters & attributes
if(0){
  View(listFilters(mart_data))
  View(listAttributes(mart_data))  
}
# Attributes are what you want returned from your mapping (i.e., what you are mapping to)
# Filters are what you are trying to map 
# Example: If you want to map Symbols to Entrez_IDs, then filters=Symbols, attributes=c(Symbols, Entrez_IDs). attributes has both because we want to return both types.
prot.map <- getBM(attributes = c("ensembl_peptide_id", "mgi_symbol"),
                  filters = "ensembl_peptide_id",
                  values = uniq.prot.bp,
                  mart = mart_data)
prot.map = prot.map[prot.map[,1] != "" & prot.map[,2] != "", ] # Remove features that didn't map
head(prot.map); nrow(prot.map)
# Check relationship of mapping 
table(table(prot.map[,1]))
table(table(prot.map[,2]))
# One-to-one mapping
# Remove Genes that map to multiple peptide IDs
prot.map = prot.map[!prot.map[,2] %in% names(table(prot.map[,2]))[table(prot.map[,2]) > 1],]
# Remove proteins from bipartite network that didn't map to a gene symbol
bp2 = bp2[bp2[,2] %in% prot.map[,1],]
# Map ensembl peptide IDs of bipartite graph to Gene Symbols
id = match(bp2[,2], prot.map[,1])
bp2[!is.na(id),2] = prot.map[id[!is.na(id)],2]

# Get g.path as edge list with edge weights
el.path = igraph::as_edgelist(g.path)
el.path = cbind(el.path, E(g.path)$weight)
# Build heterogeneous graph
el.tmp = rbind(bp2, mmi, el.path)
g.het2 = igraph::graph_from_edgelist(el = el.tmp[,1:2], directed = FALSE)
E(g.het2)$weight = as.numeric(el.tmp[,3])
# Set node type (PubChemIDs don't have any characters, and all gene symbols have at least one character)
all_chars = "[ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz]"
# sum(!grepl(all_chars, V(g.het2)$name)) == sum(!grepl(all_chars, unique(c(bp2[,1], mmi[,1], mmi[,2]))))
V(g.het2)$node_type = ifelse(!grepl(all_chars, V(g.het2)$name), "meta", "prot")

# Take largest connected component
g.het2 = largest_connected_component(g.het2)

V(g.het2)$label = paste(V(g.het2)$name, V(g.het2)$node_type, sep = "|")

#=== Multiplex-Heterogeneous and Multiplex-Homogeneous Networks: Physical+Functional PPI ----
# Physical PPI layer
phys.tmp = igraph::induced_subgraph(g.het1, which(V(g.het1)$node_type == "prot"))
phys.tmp = largest_connected_component(phys.tmp)
V(phys.tmp)$node_type = "prot_phys"
# Functional PPI layer
fun.tmp = igraph::induced_subgraph(g.het2, which(V(g.het2)$node_type == "prot"))
fun.tmp = largest_connected_component(fun.tmp)
V(fun.tmp)$node_type = "prot_fun"
# Bipartite connections
bp2[,1] = paste(bp2[,1], "meta", sep = "|")
bp2[,2] = paste(bp2[,2], "prot", sep = "|")
bp.tmp = igraph::graph_from_edgelist(el = bp2[,1:2], directed = FALSE)
E(bp.tmp)$weight = as.numeric(bp2[,3])
V(bp.tmp)$node_type = extract_string(V(bp.tmp)$name, "\\|", 2)
V(bp.tmp)$name = extract_string(V(bp.tmp)$name, "\\|", 1)
# MMI component
mmi.tmp1 = igraph::induced_subgraph(g.het1, which(V(g.het1)$node_type == "meta"))
mmi.el1 = igraph::as_edgelist(mmi.tmp1)
mmi.el1 = cbind(mmi.el1, E(mmi.tmp1)$weight)
mmi.tmp2 = igraph::induced_subgraph(g.het2, which(V(g.het2)$node_type == "meta"))
mmi.el2 = igraph::as_edgelist(mmi.tmp2)
mmi.el2 = cbind(mmi.el2, E(mmi.tmp2)$weight)
mmi.el = rbind(mmi.el1, mmi.el2)
mmi.tmp = igraph::graph_from_edgelist(as.matrix(mmi.el[,1:2]), directed = FALSE)
igraph::E(mmi.tmp)$weight = as.numeric(mmi.el[,3])
mmi.tmp = igraph::simplify(graph = mmi.tmp, edge.attr.comb = list(weight = "min"))
V(mmi.tmp)$node_type = "meta"

# Create input list
input_graph = list(prot_phys = phys.tmp,
                   prot_fun = fun.tmp,
                   "meta;prot" = bp.tmp,
                   meta = mmi.tmp)
input_graph = lapply(input_graph, function(x){
  V(x)$score = 1
  x
})
# Create multiplex-heterogeneous graph
g.mh = create_integrated_graph(graph = input_graph, data = "score", heterogeneous = TRUE, multiplex = TRUE, lcc = TRUE)

# Create new vertex attr that contains just name|component info
primary = "prot_phys"
comp = extract_string(primary, "_", 1)
tmp.nm = paste(get.type(V(g.mh)$name, 1), comp, sep = "|")
tmp.id = which(get.type(V(g.mh)$name, 3) == comp)
igraph::vertex_attr(g.mh, "agg_name", tmp.id) = tmp.nm[tmp.id]
igraph::vertex_attr(g.mh, "agg_name", -tmp.id) = V(g.mh)$name[-tmp.id]

# Multiplex-homogeneous
input_graph = list(prot_phys = phys.tmp,
                   prot_fun = fun.tmp)
input_graph = lapply(input_graph, function(x){
  V(x)$score = 1
  x
})
g.multi.homo = create_integrated_graph(graph = input_graph, data = "score", heterogeneous = FALSE, multiplex = TRUE, lcc = TRUE)

# Create new vertex attr that contains just name|component info
primary = "prot_phys"
comp = extract_string(primary, "_", 1)
tmp.nm = paste(get.type(V(g.multi.homo)$name, 1), comp, sep = "|")
tmp.id = which(get.type(V(g.multi.homo)$name, 3) == comp)
igraph::vertex_attr(g.multi.homo, "agg_name", tmp.id) = tmp.nm[tmp.id]
igraph::vertex_attr(g.multi.homo, "agg_name", -tmp.id) = V(g.multi.homo)$name[-tmp.id]

#=== Create KEGG pathway objects ===#
mh.names = extract_string(V(g.mh)$name, "\\|", 1)
mhomo.names = extract_string(V(g.multi.homo)$name, "\\|", 1)
h1.names = V(g.het1)$name
h2.names = V(g.het2)$name
kp.mh = lapply(kegg_pathways, function(x){
  # Convert Comp Name to PubChemID
  id = match(x$COMPOUND, comp.map$Query)
  tmp = comp.map$PubChem[id[!is.na(id)]]
  names(tmp) = names(x$COMPOUND)[!is.na(id)]
  x$COMPOUND = tmp
  # Remove features not in network
  x$COMPOUND = x$COMPOUND[x$COMPOUND %in% mh.names]
  x$GENE = x$GENE[x$GENE %in% mh.names]
  x
})
kp.mhomo = lapply(kegg_pathways, function(x){
  # Convert Comp Name to PubChemID
  id = match(x$COMPOUND, comp.map$Query)
  tmp = comp.map$PubChem[id[!is.na(id)]]
  names(tmp) = names(x$COMPOUND)[!is.na(id)]
  x$COMPOUND = tmp
  # Remove features not in network
  x$COMPOUND = x$COMPOUND[x$COMPOUND %in% mhomo.names]
  x$GENE = x$GENE[x$GENE %in% mhomo.names]
  x
})
kp.h1 = lapply(kegg_pathways, function(x){
  # Convert Comp Name to PubChemID
  id = match(x$COMPOUND, comp.map$Query)
  tmp = comp.map$PubChem[id[!is.na(id)]]
  names(tmp) = names(x$COMPOUND)[!is.na(id)]
  x$COMPOUND = tmp
  # Remove features not in network
  x$COMPOUND = x$COMPOUND[x$COMPOUND %in% h1.names]
  x$GENE = x$GENE[x$GENE %in% h1.names]
  x
})
kp.h2 = lapply(kegg_pathways, function(x){
  # Convert Comp Name to PubChemID
  id = match(x$COMPOUND, comp.map$Query)
  tmp = comp.map$PubChem[id[!is.na(id)]]
  names(tmp) = names(x$COMPOUND)[!is.na(id)]
  x$COMPOUND = tmp
  # Remove features not in network
  x$COMPOUND = x$COMPOUND[x$COMPOUND %in% h2.names]
  x$GENE = x$GENE[x$GENE %in% h2.names]
  x
})

#======================================#
# Construct k-folds within each pathway ----
#======================================#
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

# Access IDs of fold i for set j 
# gene_sets[[j]][fold_ids[[j]] %in% i]
# k.folds = 5
# kegg_cv = kp.mh[unlist(lapply(kp.mh, function(x) length(x$GENE) >= k.folds && length(x$COMPOUND) >= k.folds))]
# gene_fold_ids = create_folds(kfolds = k.folds, gene_sets = lapply(kegg_cv, function(x) unname(x$GENE)))
# comp_fold_ids = create_folds(kfolds = k.folds, gene_sets = lapply(kegg_cv, function(x) unname(x$COMPOUND)))
# length(gene_fold_ids)

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
tmat.res.path = paste0(path.results, "transition_matrix_multiplex_heterogeneous.rds")
rank.res.path = paste0(path.results, "node.ranking.results_multiplex_heterogeneous.rds")
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
tmat.res.path = paste0(path.results, "transition_matrix_multiplex_homogeneous.rds")
rank.res.path = paste0(path.results, "node.ranking.results_multiplex_homogeneous.rds")
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
tmat.res.path = paste0(path.results, "transition_matrix_monoplex_heterogeneous1.rds")
rank.res.path = paste0(path.results, "node.ranking.results_monoplex_heterogeneous1.rds")
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
tmat.res.path = paste0(path.results, "transition_matrix_monoplex_heterogeneous2.rds")
rank.res.path = paste0(path.results, "node.ranking.results_monoplex_heterogeneous2.rds")
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

#========#
# Figures ----
#========#
# https://wilkelab.org/cowplot/articles/drawing_with_on_plots.html 

mh.res = readRDS(paste0(path.results, "node.ranking.results_multiplex_heterogeneous.rds"))
multi.homo.res = readRDS(paste0(path.results, "node.ranking.results_multiplex_homogeneous.rds"))
h1.res = readRDS(paste0(path.results, "node.ranking.results_monoplex_heterogeneous1.rds"))
h2.res = readRDS(paste0(path.results, "node.ranking.results_monoplex_heterogeneous2.rds"))

# summary(unlist(multi.homo.res[[i]]))

# Save Dimensions
# w/o legend: w=1687,h=1282
# w/ legend: w=2500,h=1282

delta = 30
lw = 0.9
y_label = c("CDF", "Cumulative Probability")[2]

if(0){
  display.brewer.all()
  display.brewer.pal(n = 5, name = "Set1")
  brewer.pal(n = 5, name = "Set1")
}
disc.color.theme = "Set1"
div.color.theme = "RdBu"
grad.color.theme ="Blues"

#=== Multiplex-heterogeneous ----
ranks = unlist(mh.res)
groups = c()
names(mh.res) = c("RWR-MH: Low CT", "RWR-MH: Medium CT", "RWR-MH: High CT", 'RWR-MH: Mixed CT', "RWR")
for(i in seq_along(mh.res)){
  for(j in seq_along(mh.res[[i]])){
    groups = c(groups, rep(names(mh.res)[i], length(unlist(mh.res[[i]][[j]]))))
  }
}
dat = data.frame(mr = ranks, Methods = groups, row.names = NULL)

p.mh = ggplot(dat, aes(mr, colour = Methods)) +
  stat_ecdf(geom = "step", pad = FALSE, show.legend = TRUE, linewidth=lw) +
  labs(title = "Multiplex-Heterogeneous Graph",
       x = "Rank", y = y_label) +
  scale_color_manual(values = brewer.pal(n = 5, name = disc.color.theme)[c(1,2,3,5,4)]) + 
  # guides(color = "none") +
  theme(plot.title = element_text(size = 15+delta, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 13+delta, face = "bold"),
        axis.text = element_text(size = 12+delta, face = "bold.italic"),
        legend.title = element_text(size = 12+delta, face = "bold.italic"),
        legend.text = element_text(size = 11+delta, face = "italic"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major.y = element_line(colour = "black", linetype = "dashed", linewidth = 0.25))
# p.mh
if(1){
  x.lb = 0
  i=2
  x.ub = c(50,200,1000)[i]
  y.ub = c(0.18,0.4,1)[i]
  y.lb = 0
  inset = p.mh +
    labs(title=NULL) +
    coord_cartesian(xlim = c(x.lb, x.ub), ylim = c(y.lb,y.ub)) +
    theme(legend.position="none",
          axis.title = element_text(size = 4+delta, face = "bold"),
          axis.text = element_text(size = 8+delta, face = "bold.italic"),
          panel.grid.major.y = element_blank())

  ggdraw(p.mh + guides(color = "none")) +
    draw_plot(inset, x=0.35,y=0.15,width=0.5,height=0.5)
}

#=== Multiplex-homogeneous ----
ranks = unlist(multi.homo.res)
groups = c()
names(multi.homo.res) = c("RWR-MH: Low CT", "RWR-MH: Medium CT", "RWR-MH: High CT", "RWR")
for(i in seq_along(multi.homo.res)){
  for(j in seq_along(multi.homo.res[[i]])){
    groups = c(groups, rep(names(multi.homo.res)[i], length(unlist(multi.homo.res[[i]][[j]]))))
  }
}
dat = data.frame(mr = ranks, Methods = groups, row.names = NULL)

p.mho = ggplot(dat, aes(mr, colour = Methods)) +
  stat_ecdf(geom = "step", pad = FALSE, show.legend = TRUE, linewidth=lw) +
  labs(title = "Multiplex-Homogeneous Graph",
       x = "Rank", y = y_label) +
  scale_color_manual(values = brewer.pal(n = 4, name = disc.color.theme)) + 
  theme(plot.title = element_text(size = 15+delta, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 13+delta, face = "bold"),
        axis.text = element_text(size = 12+delta, face = "bold.italic"),
        legend.title = element_text(size = 12+delta, face = "bold.italic"),
        legend.text = element_text(size = 12+delta, face = "bold.italic"),
        legend.key.size = unit(1, units='cm'),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major.y = element_line(colour = "black", linetype = "dashed", linewidth = 0.25)) +
  guides(color = guide_legend(override.aes = list(linewidth = 4)))
# p.mho
if(1){
  x.lb = 0
  i=2
  x.ub = c(50,100,250,1000)[i]
  y.ub = c(0.18,0.25,0.25,1)[i]
  y.lb = 0
  inset = p.mho +
    labs(title=NULL) +
    coord_cartesian(xlim = c(x.lb, x.ub), ylim = c(y.lb,y.ub)) +
    theme(legend.position="none",
          axis.title = element_text(size = 4+delta, face = "bold"),
          axis.text = element_text(size = 8+delta, face = "bold.italic"),
          panel.grid.major.y = element_blank())
  
  ggdraw(p.mho) +
    draw_plot(inset, x=0.17,y=0.15,width=0.5,height=0.5)
}

#=== Monoplex-heterogeneous 1 ----
# Physical PPI
if(0){
  ranks = unlist(h1.res)
  groups = c()
  names(h1.res) = c("RWR-MH: Low CT", "RWR-MH: Medium CT", "RWR-MH: High CT", "RWR")
  for(i in seq_along(h1.res)){
    for(j in seq_along(h1.res[[i]])){
      groups = c(groups, rep(names(h1.res)[i], length(unlist(h1.res[[i]][[j]]))))
    }
  }
  dat = data.frame(mr = ranks, Methods = groups, row.names = NULL)
  
  p.h1 = ggplot(dat, aes(mr, colour = Methods)) +
    stat_ecdf(geom = "step", pad = FALSE, show.legend = TRUE, linewidth=lw) +
    labs(title = "Monoplex-Heterogeneous Graph",
         x = "Rank", y = y_label) +
    scale_color_manual(values = brewer.pal(n = 4, name = disc.color.theme)) + 
    theme(plot.title = element_text(size = 15+delta, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 13+delta, face = "bold"),
          axis.text = element_text(size = 12+delta, face = "bold.italic"),
          legend.title = element_text(size = 12+delta, face = "bold"),
          legend.text = element_text(size = 11+delta, face = "italic"),
          # panel.grid.major.y = element_blank(),
          panel.grid.major.y = element_line(colour = "black", linetype = "dashed", linewidth = 0.25),
          panel.background = element_rect(fill = "white", colour = "black")) +
    guides(color = "none")
  p.h1
  if(1){
    x.lb = 0
    i=3
    x.ub = c(50,100,250,1000)[i]
    y.ub = c(0.18,0.2,0.6,1)[i]
    y.lb = 0
    inset = p.h1 +
      labs(title=NULL) +
      coord_cartesian(xlim = c(x.lb, x.ub), ylim = c(y.lb,y.ub)) +
      theme(legend.position="none",
            axis.title = element_text(size = 8+delta, face = "bold"),
            axis.text = element_text(size = 8+delta, face = "bold.italic"),
            panel.grid.major.y = element_blank())
    
    ggdraw(p.h1) +
      draw_plot(inset, x=0.35,y=0.15,width=0.5,height=0.5)
  }
}

#=== Monoplex-heterogeneous 2 ----
# Functional PPI
ranks = unlist(h2.res)
groups = c()
names(h2.res) = c("RWR-MH: Low CT", "RWR-MH: Medium CT", "RWR-MH: High CT", "RWR")
for(i in seq_along(h2.res)){
  for(j in seq_along(h2.res[[i]])){
    groups = c(groups, rep(names(h2.res)[i], length(unlist(h2.res[[i]][[j]]))))
  }
}
dat = data.frame(mr = ranks, Methods = groups, row.names = NULL)

p.h2 = ggplot(dat, aes(mr, colour = Methods)) +
  stat_ecdf(geom = "step", pad = FALSE, show.legend = TRUE, linewidth=lw) +
  labs(title = "Monoplex-Heterogeneous Graph",
       x = "Rank", y = y_label) +
  scale_color_manual(values = brewer.pal(n = 4, name = disc.color.theme)) + 
  theme(plot.title = element_text(size = 15+delta, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 13+delta, face = "bold"),
        axis.text = element_text(size = 12+delta, face = "bold.italic"),
        legend.title = element_text(size = 12+delta, face = "bold"),
        legend.text = element_text(size = 11+delta, face = "italic"),
        # panel.grid.major.y = element_blank(),
        panel.grid.major.y = element_line(colour = "black", linetype = "dashed", linewidth = 0.25),
        panel.background = element_rect(fill = "white", colour = "black")) +
  guides(color = "none")
# p.h2
if(1){
  x.lb = 0
  i=3
  x.ub = c(50,100,250,1000)[i]
  y.ub = c(0.18,0.2,0.4,1)[i]
  y.lb = 0
  inset = p.h2 +
    labs(title=NULL) +
    coord_cartesian(xlim = c(x.lb, x.ub), ylim = c(y.lb,y.ub)) +
    theme(legend.position="none",
          axis.title = element_text(size = 4+delta, face = "bold"),
          axis.text = element_text(size = 8+delta, face = "bold.italic"),
          panel.grid.major.y = element_blank())
  
  ggdraw(p.h2) +
    draw_plot(inset, x=0.35,y=0.15,width=0.5,height=0.5)
}

#======================#
# RWR-MH vs. MultiXrank ----
#======================#
# RWR-MH vs. MultiXrank
if(0){
  # Apply switch probs to adj mat, then column normalize as whole?
  # Or normalize each sub-adj mat independently, then apply switch probs?
  
  # Simulate a wider set of scenarios for d, L, and aij, and discover any trends.
  # Put in Supplemental Material of paper (?)
  
  #=== Multiplex Intra/Inter layer transition probabilities ===#
  if(1){
    library(ggplot2)
    library(cowplot)
    library(grid)
    library(RColorBrewer)
    
    # Save Dimensions
    # Intra: w=1687,h=1282
    # Inter: w=2000,h=1282
    
    delta = 30
    lw = 1.4
    
    a = seq(0, 1, 0.005)
    # node degree
    d = 50 # The sensitivity of trans probs to switch prob decreases as degree of node increases.
    # number of layers
    L = 3
    # edge weight
    aij = 0.8
    # MultiXrank: exponential
    mxr.intra = aij / (d + a/(1-a))
    mxr.inter = (1/(L-1)) / (d*(1-a)/a + 1)
    # RWR-MH: linear
    rwr.intra = (aij / d) * (1-a)
    rwr.inter = a / (L-1)
    X = matrix(rep(a, 2), ncol = 2)
    Y1 = matrix(c(mxr.intra, rwr.intra), ncol = 2)
    Y2 = matrix(c(mxr.inter, rwr.inter), ncol = 2)
    
    dat = data.frame(X = rep(a, 2),
                     intra = c(mxr.intra, rwr.intra),
                     inter = c(mxr.inter, rwr.inter),
                     group = c(rep("mxr",length(mxr.inter)), rep("rwrmh",length(rwr.inter))))
    unique(sort(dat$group))
    
    p.intra = ggplot(dat, aes(x=X, y=intra)) +
      geom_line(aes(color=group), linewidth=lw) +
      scale_color_manual(values = brewer.pal(n = 3, name = disc.color.theme)[1:2],
                         labels = c("MultiXrank", "RWR-MH")) +
      labs(title = "Multiplex Intra-Layer Adjacency Matrix",
           x = "Switch Probability",
           y = "Transition Probability",
           color = "") +
      theme_light() +
      theme(plot.title = element_text(size=14+delta,face='bold', hjust=0.5),
            axis.text = element_text(size=12+delta,face='bold.italic'),
            axis.title.x = element_text(size=12+delta,face='bold',margin=margin(t=15)),
            axis.title.y = element_text(size=12+delta,face='bold',margin=margin(r=15)),
            panel.grid.major.y=element_line(colour = "black", linetype = "dashed", linewidth = 0.25),
            panel.grid.major.x=element_blank(),
            panel.grid.minor=element_blank(),
            legend.text = element_text(size=14+delta,face='bold.italic')) +
      guides(color = "none")
    
    p.inter = ggplot(dat, aes(x=X, y=inter)) +
      geom_line(aes(color=group), linewidth=lw) +
      scale_color_manual(values = brewer.pal(n = 3, name = disc.color.theme)[1:2],
                         labels = c("MultiXrank", "RWR-MH")) +
      labs(title = "Multiplex Inter-Layer Adjacency Matrix",
           x = "Switch Probability",
           y = "Transition Probability",
           color = "") +
      theme_light() +
      theme(plot.title = element_text(size=14+delta,face='bold', hjust=0.5),
            axis.text = element_text(size=12+delta,face='bold.italic'),
            axis.title.x = element_text(size=12+delta,face='bold',margin=margin(t=15)),
            axis.title.y = element_text(size=12+delta,face='bold',margin=margin(r=15)),
            panel.grid.major.y=element_line(colour = "black", linetype = "dashed", linewidth = 0.25),
            panel.grid.major.x=element_blank(),
            panel.grid.minor=element_blank(),
            legend.text = element_text(size=14+delta,face='bold.italic')) +
      guides(color = guide_legend(override.aes = list(linewidth = 4)))
    
    #=== Bipartite adj mat (transition from component alpha to beta) ===#
    mxr = TRUE
    if(mxr){ # MultiXrank assumptions
      bij = 0.8
      bp.deg_alpha.j = 10 # number of bipartite neighbors of node j, component alpha (any layer) that are in component beta (any layer)
      L.beta = 4 # number of layers in component beta
      psi_j = 2 # number of other components that node j of current layer of component alpha is connected to.
      L.beta_alpha.j = L.beta - 0 # number of layers in component beta that node j of current layer of component alpha is connected to.
      mxr.bp = a * bij / (bp.deg_alpha.j * L.beta)
      rwr.bp = a * bij / (bp.deg_alpha.j * psi_j * L.beta_alpha.j)
    }else{ # RWR-MH assumptions
      bij = 0.8
      bp.deg_alpha.j.mxr = 10 # number of bipartite neighbors of node j, component alpha (any layer) that are in component beta (any layer)
      L.beta = 4 # number of layers in component beta
      bp.deg_alpha.j.rwr = 5 # number of bipartite neighbors of node j, component alpha (current layer) that are in component beta (specific layer)
      psi_j = 2 # number of other components that node j of current layer of component alpha is connected to.
      L.beta_alpha.j = L.beta / 2 # number of layers in component beta that node j of current layer of component alpha is connected to.
      mxr.bp = a * bij / (bp.deg_alpha.j.mxr * L.beta)
      rwr.bp = a * bij / (bp.deg_alpha.j.rwr * psi_j * L.beta_alpha.j)
    }
    dat = data.frame(X = rep(a, 2),
                     Y = c(mxr.bp, rwr.bp),
                     group = c(rep("mxr", length(mxr.bp)), rep("rwrmh", length(rwr.bp))))
    
    p.bp = ggplot(dat, aes(x=X, y=Y)) +
      geom_line(aes(color=group), linewidth=1.2) +
      scale_color_manual(values = brewer.pal(n = 3, name = disc.color.theme)[1:2],
                         labels = c("MultiXrank", "RWR-MH")) +
      labs(title = "Multiplex Inter-Layer Adjacency Matrix",
           x = "Jump Probability",
           y = "Transition Probability",
           color = "") +
      theme_light() +
      theme(plot.title = element_text(face='bold'),
            axis.text = element_text(size=12,face='bold.italic'),
            axis.title.x = element_text(face='bold',margin=margin(t=15)),
            axis.title.y = element_text(face='bold',margin=margin(r=15)),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            legend.text = element_text(size=12))
    
    matplot(x = X, y = Y, type = "b", pch = 1, col = 1:2, main = "Bipartite Adjacency Matrix", xlab = "Jump Probability", ylab = "Transition Probability")
    legend("bottomright", legend = c("MultiXrank", "RWR-MH"), col = 1:2, pch = 1)
    # As L.beta increases, MultiXrank and RWR-MH bipartite trans probs become less sensitive to jump prob.
    # As psi_j increases, RWR-MH bipartite trans probs become less sensitive to jump prob.
  }
  #=====END=====#
  # Old code
  a = seq(0.001, 0.999, 0.005)
  if(0){ # Multiplex Intra/Inter
    d = 50 # The sensitivity of trans probs to switch prob decreases as degree of node increases.
    L = 3
    aij = 0.8
    # MultiXrank: exponential
    mxr.intra = aij / (d + a/(1-a))
    mxr.inter = (1/(L-1)) / (d*(1-a)/a + 1)
    # RWR-MH: linear
    rwr.intra = (aij / d) * (1-a)
    rwr.inter = a / (L-1)
    X = matrix(rep(a, 2), ncol = 2)
    Y1 = matrix(c(mxr.intra, rwr.intra), ncol = 2)
    Y2 = matrix(c(mxr.inter, rwr.inter), ncol = 2)
    par(mfrow = c(1,2))
    matplot(x = X, y = Y1, type = "b", pch = 1, col = 1:2, main = "Multiplex within-layer Adjacency Matrix", xlab = "Switch Probability", ylab = "Transition Probability")
    legend("bottomleft", legend = c("MultiXrank", "RWR-MH"), col = 1:2, pch = 1)
    matplot(x = X, y = Y2, type = "b", pch = 1, col = 1:2, main = "Multiplex between-layer Adjacency Matrix", xlab = "Switch Probability", ylab = "Transition Probability")
    legend("topleft", legend = c("MultiXrank", "RWR-MH"), col = 1:2, pch = 1)
  }
  if(0){ # Bipartite adj mat (transition from component alpha to beta)
    par(mfrow = c(1,1))
    mxr = TRUE
    if(mxr){ # MultiXrank assumptions
      bij = 0.8
      bp.deg_alpha.j = 10 # number of bipartite neighbors of node j, component alpha (any layer) that are in component beta (any layer)
      L.beta = 4 # number of layers in component beta
      psi_j = 2 # number of other components that node j of current layer of component alpha is connected to.
      L.beta_alpha.j = L.beta - 0 # number of layers in component beta that node j of current layer of component alpha is connected to.
      mxr.bp = a * bij / (bp.deg_alpha.j * L.beta)
      rwr.bp = a * bij / (bp.deg_alpha.j * psi_j * L.beta_alpha.j)
    }else{ # RWR-MH assumptions
      bij = 0.8
      bp.deg_alpha.j.mxr = 10 # number of bipartite neighbors of node j, component alpha (any layer) that are in component beta (any layer)
      L.beta = 4 # number of layers in component beta
      bp.deg_alpha.j.rwr = 5 # number of bipartite neighbors of node j, component alpha (current layer) that are in component beta (specific layer)
      psi_j = 2 # number of other components that node j of current layer of component alpha is connected to.
      L.beta_alpha.j = L.beta / 2 # number of layers in component beta that node j of current layer of component alpha is connected to.
      mxr.bp = a * bij / (bp.deg_alpha.j.mxr * L.beta)
      rwr.bp = a * bij / (bp.deg_alpha.j.rwr * psi_j * L.beta_alpha.j)
    }
    
    X = matrix(rep(a, 2), ncol = 2)
    Y = matrix(c(mxr.bp, rwr.bp), ncol = 2)
    matplot(x = X, y = Y, type = "b", pch = 1, col = 1:2, main = "Bipartite Adjacency Matrix", xlab = "Jump Probability", ylab = "Transition Probability")
    legend("bottomright", legend = c("MultiXrank", "RWR-MH"), col = 1:2, pch = 1)
    # As L.beta increases, MultiXrank and RWR-MH bipartite trans probs become less sensitive to jump prob.
    # As psi_j increases, RWR-MH bipartite trans probs become less sensitive to jump prob.
  }
}




