#==========================================#
# Data for Gene Prioritization Task:
#   OMIM, GO (cc, bp, mf), and STRING PPIN 
#==========================================#
# Attach necessary libraries
library(mgsa)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(igraph)
library(biomaRt)

# Path names
path.data = "/path/to/data/"
path.results = "/path/to/results/"

#========================#
# PPI Network from STRING ----
#========================#
# Proteins as Ensembl peptide IDs 
edge.weight.cutoff = 0.5

ppi_el = data.table::fread(file = paste0(path.data, "9606.protein.links.v12.0.txt.gz"),
                           header = TRUE, sep = " ") %>% 
  dplyr::mutate(combined_score = combined_score / 1000) %>%
  dplyr::filter(combined_score >= edge.weight.cutoff) %>%
  dplyr::mutate(protein1 = do.call("c", lapply(strsplit(protein1, "\\."), function(x) x[2])),
                protein2 = do.call("c", lapply(strsplit(protein2, "\\."), function(x) x[2]))) %>%
  dplyr::distinct() %>%
  as.data.frame()

# Remove duplicate edges
ppi_mat = as.matrix(ppi_el[,1:2])
temp = character(nrow(ppi_mat))
for(i in 1:length(temp)){ 
  if(i %% 100000 == 0) message(i)
  x = sort(c(ppi_mat[i,1], ppi_mat[i,2]))
  temp[i] = paste(x, collapse = "|")
}

ppi_el = data.frame(proteins = temp, combined_score = ppi_el$combined_score) %>%
  dplyr::distinct() %>%
  dplyr::mutate(protein1 = do.call("c", lapply(strsplit(proteins, "\\|"), function(x) x[1])),
                protein2 = do.call("c", lapply(strsplit(proteins, "\\|"), function(x) x[2]))) %>%
  dplyr::select(protein1, protein2, combined_score)

g = igraph::graph_from_edgelist(el = as.matrix(ppi_el[,1:2]), directed = FALSE)
E(g)$weight = ppi_el$combined_score

if(!igraph::is_connected(g)){
  comps = igraph::components(g)
  largest_comp_id = which.max(comps$csize)
  g = igraph::induced_subgraph(g, which(comps$membership == largest_comp_id)) 
}

#=================#
# Mapping Gene IDs ----
#=================#
# listMarts()
mart.obj = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://jul2023.archive.ensembl.org")
# listDatasets(mart = mart.obj)
mart.obj = useDataset(dataset = "hsapiens_gene_ensembl", mart = mart.obj)
# View(listAttributes(mart = mart.obj))
# View(listFilters(mart = mart.obj))

# Mapping Ensembl peptide IDs from PPIN to Entrez IDs
map.values.ensp = V(g)$name
map.ensp = getBM(attributes = c("entrezgene_id", "ensembl_peptide_id"),
                 filters = "ensembl_peptide_id",
                 values = map.values.ensp,
                 mart = mart.obj) %>%
  dplyr::distinct() %>%
  dplyr::filter(ensembl_peptide_id != "",
                entrezgene_id != "") %>%
  as.data.frame()

# nrow(map.ensp)
# table(table(map.ensp$ensembl_peptide_id))
# sum(!map.values.ensp %in% map.ensp$ensembl_peptide_id) # 551 peptides in PPIN don't map to Entrez ID

# Removing peptides from PPIN that didn't map to an Entrez ID. This may cause the graph to become disconnected.
no.map = map.values.ensp[!map.values.ensp %in% map.ensp$ensembl_peptide_id]
g = igraph::delete_vertices(g, which(V(g)$name %in% no.map))
if(!igraph::is_connected(g)){ # Get largest connected component
  comps = igraph::components(g)
  largest_comp_id = which.max(comps$csize)
  g = igraph::induced_subgraph(g, which(comps$membership == largest_comp_id)) 
}

# Taking first mapping for every Ensembl peptide ID, giving priority to those that mapped to an Entrez ID in OMIM data
multi.map = names(table(map.ensp$ensembl_peptide_id))[table(map.ensp$ensembl_peptide_id) > 1]
rm.id = c()
for(i in seq_along(multi.map)){
  tmp.ids = which(map.ensp$ensembl_peptide_id == multi.map[i])
  cond = map.ensp$entrezgene_id[tmp.ids] %in% omim_genes_all
  if(all(!cond)){
    rm.id = c(rm.id, tmp.ids[-1])
  }else{
    rm.id = c(rm.id, tmp.ids[cond][-1], tmp.ids[!cond])
  }
}
map.ensp = map.ensp[-rm.id,]
# table(table(map.ensp$ensembl_peptide_id)) # Verify

V(g)$entrez_id = map.ensp$entrezgene_id[match(V(g)$name, map.ensp$ensembl_peptide_id)]
V(g)$ensp_id = V(g)$name
V(g)$name = V(g)$entrez_id
g = igraph::delete_vertex_attr(g, "entrez_id")

# Keep node with largest degree from duplicate Entrez IDs (total of 33 IDs)
tmp = table(V(g)$name)
ttmp = table(tmp)[as.numeric(names(table(tmp))) > 1]
N.tmp = sum((as.numeric(names(ttmp)) - 1) * ttmp)
multi.nodes = names(tmp)[tmp > 1]
all.d = igraph::degree(g)
FUN = max # Function to determine which node to KEEP.
rm.id = numeric(N.tmp)
idx = 1
for(i in 1:length(multi.nodes)){
  d = igraph::degree(g, which(V(g)$name == multi.nodes[i]))
  n.tmp = length(d)
  if(length(unique(d)) > 1){
    tmp.id = which(V(g)$name == multi.nodes[i] & all.d != FUN(d))
  }else{
    tmp.id = sample(x = which(V(g)$name == multi.nodes[i]), size = length(d) - 1)
  }
  rm.id[idx:(idx + n.tmp - 2)] = tmp.id
  idx = idx + n.tmp - 1
}
g = igraph::delete_vertices(graph = g, v = rm.id)

# Get largest connected component
if(!igraph::is_connected(g)){ 
  comps = igraph::components(g)
  largest_comp_id = which.max(comps$csize)
  g = igraph::induced_subgraph(g, which(comps$membership == largest_comp_id)) 
}

saveRDS(g, file = paste0(path.results, "STRING_0.5_igraph.rds"))
