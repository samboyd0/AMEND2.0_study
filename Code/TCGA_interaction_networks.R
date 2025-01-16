#===============================#
# Molecular Interaction Networks
#     for TCGA Analysis
#===============================#
# AMEND vignette Question: Should I cut out data creation steps if possible. Focus more on using the package?
#   Or, create a separate document with clearer examples?

rm(list = ls())

library(igraph)
library(readxl)
library(dplyr)
library(biomaRt)
library(graphite)

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
get_random_names = function(n, rng.seed=NULL){
  if(!is.null(rng.seed) && is.numeric(rng.seed)) set.seed(rng.seed)
  grid = expand.grid(letters, 0:9, toupper(letters), stringsAsFactors = FALSE)
  random.name = character(n)
  avail.ids = 1:nrow(grid)
  for(i in seq_along(random.name)){
    ids = sample(avail.ids, 2)
    avail.ids = setdiff(avail.ids, ids)
    tmp = c(as.character(grid[ids[1],]), as.character(grid[ids[2],]))
    random.name[i] = paste(tmp[sample(1:length(tmp),length(tmp))],collapse="")
  }
  random.name
}

# PPI physical interactions from STRING
# Kidney-specific PPI network
# Read papers using multiplex on PPINs with PPIs from multiple databases (Menche's paper?)
#   (Should I use another db other than STRING?)
# Aggregate multiplex layers

# Use the 'graphite' package to get functional interactions (i.e., pathway networks from annotation databases)
# Use the 'aracne.networks' package to get a GRN from TCGA-KIRC gene expression data

#================#
# Script Settings
#================#
set.seed(2707083)

ppi_edge_threshold = 0.7

path.data = "/path/to/data/"
path.amend = "/path/to/AMEND_code/"
path.amend.res = "/path/to/results/"

#==================#
# Load in TCGA Data
#==================#
# hr = readRDS("hazard_ratios.rds")
# omics.dat = readRDS("multi_omics.rds")
# de.res = readRDS("DE_results.rds")
source("TCGA_KIRC_Preprocessing.R")

# Remove everything after '|' for miRNA names
rownames(de.res$mirna) = extract_string(rownames(de.res$mirna), "\\|", 1)
names(hr$mirna) = extract_string(names(hr$mirna), "\\|", 1)

# Changing 'MIRXX-1AS' to 'MIRXXB1' for genes in methylation dataset
id = which(grepl('AS',rownames(de.res$methyl)) & grepl('MIR',rownames(de.res$methyl)))
for(i in id){
  tips = extract_string(rownames(de.res$methyl)[i],'-',1)
  ends = extract_string(rownames(de.res$methyl)[i],'-',2)
  tips = paste0(tips, "B")
  ends = gsub('AS', '', ends)
  rownames(de.res$methyl)[i] = paste0(tips, ends)
}
id = which(grepl('AS',names(hr$methyl)) & grepl('MIR',names(hr$methyl)))
for(i in id){
  tips = extract_string(names(hr$methyl)[i],'-',1)
  ends = extract_string(names(hr$methyl)[i],'-',2)
  tips = paste0(tips, "B")
  ends = gsub('AS', '', ends)
  names(hr$methyl)[i] = paste0(tips, ends)
}

#===================#
# miRNA-mRNA Network ----
#===================#
dat = read_xlsx(path = paste0(path.data, "miRTarBase_SE_WR.xlsx")) %>%
  dplyr::filter(substr(miRNA,1,3) == "hsa" & `Species (miRNA)` == "Homo sapiens" & `Species (Target Gene)` == "Homo sapiens") %>%
  dplyr::rename(ID = `miRTarBase ID`,
                Symbol = `Target Gene`,
                EntrezID = `Target Gene (Entrez ID)`) %>%
  dplyr::select(ID, miRNA, Symbol, EntrezID) %>%
  dplyr::distinct()

mirna.el = as.matrix(dat[,c('miRNA', 'EntrezID')])
# Removing leading whitespace
mirna.el[,2] = gsub(pattern=" ", replacement="", x=mirna.el[,2])
bp.mirna.mrna = igraph::graph_from_edgelist(el = mirna.el, directed = FALSE)
V(bp.mirna.mrna)$node_type = ifelse(V(bp.mirna.mrna)$name %in% mirna.el[,1], "mirna", "mrna")
E(bp.mirna.mrna)$weight = ppi_edge_threshold

vertex_attr(bp.mirna.mrna, "native_nom") = ifelse(V(bp.mirna.mrna)$node_type == "mirna", "miRID", "EntrezID")
id = which(V(bp.mirna.mrna)$node_type == "mirna")
vertex_attr(bp.mirna.mrna, "miRID", id) = V(bp.mirna.mrna)$name[id]
id = which(V(bp.mirna.mrna)$node_type == "mrna")
vertex_attr(bp.mirna.mrna, "EntrezID", id) = V(bp.mirna.mrna)$name[id]

#=============#
# PPI Networks ----
#=============#
# Species: Homo Sapiens
# 4 layers:
#   - Kidney-specific PPIs (edgelist)(Entrez)
#   - Reactome pathway network from graphite (igraph)(UniProt)
#   - GRN from ARACNE on TCGA-KIRC data (edgelist)(Entrez)
#   - Physical PPI from STRING (edgelist)(Ensembl)

#=== Kidney-specific PPIN ===#
# Source: https://snap.stanford.edu/biodata/datasets/10013/10013-PPT-Ohmnet.html 
k.ppi = data.table::fread(file = paste0(path.data, "PPT-Ohmnet_tissues-combined.txt")) %>%
  dplyr::filter(tissue == "kidney") %>%
  dplyr::select(-tissue) %>%
  as.matrix()
k.ppi = matrix(as.character(k.ppi), ncol=2, dimnames = dimnames(k.ppi))
k.ppi = cbind(k.ppi, as.character(ppi_edge_threshold))
g.ohmnet = igraph::graph_from_edgelist(el = k.ppi[,1:2], directed = FALSE)
E(g.ohmnet)$weight = as.numeric(k.ppi[,3])

vertex_attr(g.ohmnet, "native_nom") = "EntrezID"
vertex_attr(g.ohmnet, "EntrezID") = V(g.ohmnet)$name

#=== Pathway Network (graphite) ===#
# graphite::pathwayDatabases()
hsR = graphite::pathways("hsapiens", "reactome")
print(length(hsR))
for(i in seq_along(hsR)){
  if(i %% 100 == 0) message(i)
  g.tmp = graphite::pathwayGraph(pathway = hsR[[i]], which = "proteins")
  g.tmp = igraph::graph_from_graphnel(graphNEL = g.tmp)
  g.tmp = igraph::as.undirected(graph = g.tmp, mode = "collapse")
  if(i == 1) g.path = g.tmp else g.path = igraph::union(g.path, g.tmp)
}
edge.weight.attrs = igraph::edge_attr_names(g.path)[grepl("weight", igraph::edge_attr_names(g.path))]
for(i in seq_along(edge.weight.attrs)){
  g.path = igraph::delete_edge_attr(graph = g.path, name = edge.weight.attrs[i])
}
E(g.path)$weight = ppi_edge_threshold # Give all edges the minimum edge weight
# Remove nodes that don't contain 'UNIPROT:'. These will be CHEBI (chemical entities of biological interest)
g.path = igraph::delete_vertices(graph = g.path, v = which(!grepl('UNIPROT:', V(g.path)$name)))
V(g.path)$name = gsub("UNIPROT:", "", V(g.path)$name)
# saveRDS(g.path, file = "TCGA_graphite_network.rds")
# g.path = readRDS(paste0(path.data, "TCGA_graphite_network.rds"))
vertex_attr(g.path, "native_nom") = "UniProt"
vertex_attr(g.path, "UniProt") = V(g.path)$name

#=== GRN (ARACNE) ===#
# Writing ARACNE regulon object as edgelist
if(0){
  library(aracne.networks)
  data(regulonkirc)
  write.regulon(regulon = regulonkirc, 
                n = Inf,
                file = paste0(path.data, "aracne_el.txt"),
                sep = "\t", header = TRUE)
}
# Reading in ARACNE network as edgelist
grn = data.table::fread(file = paste0(path.data, "aracne_el.txt"), sep="\t", header=TRUE)
grn = grn[grn$likelihood >= ppi_edge_threshold,-3]
grn = as.matrix(grn)
grn = matrix(as.character(grn), ncol=3, dimnames = dimnames(grn))
g.aracne = igraph::graph_from_edgelist(el = grn[,1:2], directed = FALSE)
E(g.aracne)$weight = as.numeric(grn[,3])

vertex_attr(g.aracne, "native_nom") = "EntrezID"
vertex_attr(g.aracne, "EntrezID") = V(g.aracne)$name

#=== STRING Physical PPI Network ===#
ppi = data.table::fread(file = paste0(path.data, "9606.protein.physical.links.v12.0.txt.gz"), header = T) %>%
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
string.ppi = matrix(c(extract_string(dat[,1], "\\|", 1), extract_string(dat[,1], "\\|", 2), dat[,2]), ncol = 3)
g.string = igraph::graph_from_edgelist(el = string.ppi[,1:2], directed = FALSE)
E(g.string)$weight = as.numeric(string.ppi[,3])

vertex_attr(g.string, "native_nom") = "ENSPID"
vertex_attr(g.string, "ENSPID") = V(g.string)$name

#===============================#
#=== PPIN Annotation Mapping ===#
#===============================#
# list the databases available from biomaRt and archived versions
if(0){
  listMarts() 
  listEnsemblArchives()
} 
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://jan2024.archive.ensembl.org")
# List datasets available within a certain database
if(0) listDatasets(mart) 
mart_data <- useDataset("hsapiens_gene_ensembl", mart = mart)
# Check the names of available filters & attributes
if(0){
  View(listFilters(mart_data))
  View(listAttributes(mart_data))  
}
# Attributes are what you want returned from your mapping (i.e., what you are mapping to)
# Filters are what you are trying to map 

## Ensembl peptide IDs to Entrez IDs
uniq.ensp = unique(c(string.ppi[,1], string.ppi[,2]))
ensp2entrez <- getBM(attributes = c("ensembl_peptide_id", "entrezgene_id"),
                 filters = "ensembl_peptide_id",
                 values = uniq.ensp,
                 mart = mart_data)
# Remove features that didn't map
ensp2entrez = ensp2entrez[apply(ensp2entrez,1,function(x) all(x != "")), ] 
ensp2entrez[,2] = as.character(ensp2entrez[,2])

# head(ensp2entrez)
# table(table(ensp2entrez[,1]))
# table(table(ensp2entrez[,2]))
# 
# mean(uniq.ensp %in% ensp2entrez[,1]) # 0.98

## UniProt IDs to Entrez IDs
uniq.uniprot = unique(V(g.path)$name)
# 'uniprotswissprot' filter has better mapping... I used 'uniprot_gn_id' for the mmus graphite graph...
biomart.filter = c("uniprot_gn_id", "uniprotswissprot")[2]
uniprot2entrez <- getBM(attributes = c(biomart.filter, "entrezgene_id"),
                     filters = biomart.filter,
                     values = uniq.uniprot,
                     mart = mart_data)
# Remove features that didn't map
uniprot2entrez = uniprot2entrez[apply(uniprot2entrez,1,function(x) all(x != "")), ] 
uniprot2entrez[,2] = as.character(uniprot2entrez[,2])

# nrow(uniprot2entrez)
# head(uniprot2entrez)
# table(table(uniprot2entrez[,1]))
# table(table(uniprot2entrez[,2]))
# 
# mean(uniq.uniprot %in% uniprot2entrez[,1]) # 0.95

#=== Apply mappings ===#
# ensp2entrez = rbind(ensp2entrez, "0")
id = match(V(g.string)$name, ensp2entrez[,1])
V(g.string)$name[!is.na(id)] = make.unique(ensp2entrez[id[!is.na(id)],2], sep=".")
V(g.string)$EntrezID[!is.na(id)] = ensp2entrez[id[!is.na(id)],2]
# V(g.string)$name[is.na(id)] = get_random_names(sum(is.na(id)))

id = match(V(g.path)$name, uniprot2entrez[,1])
V(g.path)$name[!is.na(id)] = make.unique(uniprot2entrez[id[!is.na(id)],2], sep=".")
V(g.path)$EntrezID[!is.na(id)] = uniprot2entrez[id[!is.na(id)],2]
# V(g.path)$name[is.na(id)] = get_random_names(sum(is.na(id)))

#====================#
# Methylation Network ----
#====================#
# Map symbols of methylation and mRNA-seq data to Entrez IDs
uniq.symbols = unique(c(rownames(de.res$mrna), rownames(de.res$methyl)))

biomart.filter = c("hgnc_symbol")[1]
symbol2entrez <- getBM(attributes = c(biomart.filter, "entrezgene_id"),
                        filters = biomart.filter,
                        values = uniq.symbols,
                        mart = mart_data)
# Remove features that didn't map
symbol2entrez = symbol2entrez[apply(symbol2entrez,1,function(x) all(x != "")), ] 

# nrow(symbol2entrez)
# head(symbol2entrez)
# 
# table(table(symbol2entrez[,1]))
# table(table(symbol2entrez[,2]))
# 
# mean(rownames(de.res$mrna) %in% symbol2entrez[,1]) # 0.85
# mean(rownames(de.res$methyl) %in% symbol2entrez[,1]) # 0.82
# mean(rownames(de.res$methyl)[substr(rownames(de.res$methyl),1,3) == "MIR"] %in% symbol2entrez[,1]) # 0.92

# Add column 'EntrezID' to de.res$mrna
id = match(rownames(de.res$mrna), symbol2entrez[,1])
de.res$mrna$EntrezID[!is.na(id)] = symbol2entrez[id[!is.na(id)],2]

# Create Bipartite graph object between methylation component and PPIN
uniq.ppi = unique(c(V(g.string)$name, V(g.aracne)$name, V(g.ohmnet)$name, V(g.path)$name))

id = match(rownames(de.res$methyl), symbol2entrez[,1])
de.res$methyl$EntrezID[!is.na(id)] = symbol2entrez[id[!is.na(id)],2]

methyl.entrez.dat = de.res$methyl[!is.na(de.res$methyl$EntrezID),]
methyl.entrez.id = methyl.entrez.dat$EntrezID

common = intersect(methyl.entrez.id, uniq.ppi)

# mean(methyl.entrez.id %in% common) # 0.90
# mean(uniq.ppi %in% common) # 0.77

bp.el = matrix(c(paste(rownames(methyl.entrez.dat)[match(common, methyl.entrez.dat$EntrezID)], "methyl", sep="|"), paste(common, "mrna", sep="|")), ncol=2)
bp.methyl.mrna = igraph::graph_from_edgelist(el = bp.el, directed = FALSE)
V(bp.methyl.mrna)$node_type = extract_string(V(bp.methyl.mrna)$name,"\\|",2)
V(bp.methyl.mrna)$name = extract_string(V(bp.methyl.mrna)$name,"\\|",1)

vertex_attr(bp.methyl.mrna, "native_nom") = ifelse(V(bp.methyl.mrna)$node_type == "mrna", "EntrezID", "Symbol")
id = which(V(bp.methyl.mrna)$node_type == "mrna")
vertex_attr(bp.methyl.mrna, "EntrezID", id) = V(bp.methyl.mrna)$name[id]
id = which(V(bp.methyl.mrna)$node_type == "methyl")
vertex_attr(bp.methyl.mrna, "Symbol", id) = V(bp.methyl.mrna)$name[id]
id = which(V(bp.methyl.mrna)$node_type == "methyl" & V(bp.methyl.mrna)$name %in% rownames(methyl.entrez.dat))
vertex_attr(bp.methyl.mrna, 'EntrezID', id) = methyl.entrez.dat$EntrezID[match(V(bp.methyl.mrna)$name[id], rownames(methyl.entrez.dat))]

# miRNA nomenclature resource: 
# https://old.abmgood.com/marketing/knowledge_base/miRNA_Introduction.php#:~:text=Instead%20of%20the%20miR%2FmiR,original%20miR%2FmiR*%20names. 

#=== miRNA names from Methylation data ===#
# miRNA from methylation data
id=grep("MIR", rownames(de.res$methyl))
mirna.meth.nm = rownames(de.res$methyl)[id]
if(0){
  mirna.meth.nm[grep("-",mirna.meth.nm)]
  mirna.meth.nm[!grepl("-",mirna.meth.nm)]
  unique(substr(mirna.meth.nm,1,3))
  any(mirna.meth.nm[!grepl("-",mirna.meth.nm)] %in% extract_string(mirna.meth.nm[grep("-",mirna.meth.nm)],"-",1))
  mirna.meth.nm[grepl("let",mirna.meth.nm,ignore.case = T)]
}

#=== miRNA-seq Names ===#
# Currently only mature miRNA ('miR' rather than 'mir')
mirna.de.nm = rownames(de.res$mirna)
if(0){
  unique(substr(mirna.de.nm,1,7))
  ends = substr(mirna.de.nm, 9, 1000)
  end.len = unlist(lapply(strsplit(ends, "-"), length))
  table(end.len)
  ends[end.len == 3]
  unique(extract_string(ends[end.len == 2],"-",2))
  mirna.de.nm[grepl('let',mirna.de.nm,ignore.case = T)]
  
  # Remove 5p and 3p
  gsub(pattern='-[35]p',replacement='',x=ends)
  
  # Check for 'as'
  mirna.de.nm[grepl('as',mirna.de.nm,ignore.case = T)]
}
if(0){
  ### Create aliases of miRNA-seq names that will accord with miRNA gene names ###
  ## Drop '5p' & '3p'
  new.names = gsub(pattern = '-[35]p',replacement = '', x = mirna.de.nm)
  ## Replace 'hsa-miR' with 'MIR' and 'hsa-let' with 'MIRLET'
  new.names = gsub(pattern = 'hsa-miR', replacement = 'MIR', x = new.names)
  new.names = gsub(pattern = 'hsa-let', replacement = 'MIRLET', x = new.names)
  ## Capitalize all letters
  new.names = toupper(new.names)
  ## Remove last hyphen from 'MIR-XXA-1' ('XX' followed by letter followed by hyphen)
  tips = extract_string(new.names,"-",1)
  ends = unlist(lapply(strsplit(new.names, "-"), function(x) paste(x[2:length(x)],collapse="-") ))
  hyphen = grepl('-',ends)
  ABC.b4.hyphen = grepl('[ABCDEFGHIJKLMNOPQRSTUVWXYZ]', extract_string(ends,'-',1))
  ends[hyphen & ABC.b4.hyphen] = gsub(pattern = '-', replacement = '', x = ends[hyphen & ABC.b4.hyphen])
  new.names = paste0(tips, ends)
}

#=== miRNA BP graph Names ===#
# Currently only mature miRNA ('miR' rather than 'mir')
mirna.g.nm = V(bp.mirna.mrna)$name[V(bp.mirna.mrna)$node_type == "mirna"]
if(0){
  unique(substr(mirna.g.nm,1,7))
  ends = substr(mirna.g.nm, 9, 1000)
  end.len = unlist(lapply(strsplit(ends, "-"), length))
  table(end.len)
  ends[end.len == 3]
  unique(extract_string(ends[end.len == 2],"-",2))
  mirna.g.nm[grepl('let',mirna.g.nm,ignore.case = T)]
  
  # Remove 5p and 3p
  gsub(pattern='-[35]p',replacement='',x=ends)
  
  # Check for 'as'
  mirna.g.nm[grepl('as',mirna.g.nm,ignore.case = T)]
}

### Create aliases of miRNA-seq names that will accord with miRNA gene names ###
## Drop '5p' & '3p'
new.names = gsub(pattern = '-[35]p',replacement = '', x = mirna.g.nm)
## Replace 'hsa-miR' with 'MIR' and 'hsa-let' with 'MIRLET'
new.names = gsub(pattern = 'hsa-miR', replacement = 'MIR', x = new.names)
new.names = gsub(pattern = 'hsa-let', replacement = 'MIRLET', x = new.names)
## Capitalize all letters
new.names = toupper(new.names)
## Remove last hyphen from 'MIR-XXA-1' ('XX' followed by letter followed by hyphen)
tips = extract_string(new.names,"-",1)
ends = unlist(lapply(strsplit(new.names, "-"), function(x) paste(x[2:length(x)],collapse="-") ))
hyphen = grepl('-',ends)
ABC.b4.hyphen = grepl('[ABCDEFGHIJKLMNOPQRSTUVWXYZ]', extract_string(ends,'-',1))
ends[hyphen & ABC.b4.hyphen] = gsub(pattern = '-', replacement = '', x = ends[hyphen & ABC.b4.hyphen])
new.names = paste0(tips, ends)

common = intersect(new.names, rownames(de.res$methyl))

# mean(new.names %in% common)
# mean(mirna.meth.nm %in% common)

bp.el = NULL
for(i in seq_along(common)){
  N = sum(new.names == common[i])
  tmp.nm = mirna.g.nm[new.names %in% common[i]]
  tmp = matrix(c(paste(tmp.nm,'mirna',sep='|'), paste(rep(common[i], N),'methyl',sep='|')),ncol=2)
  bp.el = rbind(bp.el, tmp)
}

bp.methyl.mirna = igraph::graph_from_edgelist(el = bp.el, directed = FALSE)
V(bp.methyl.mirna)$node_type = extract_string(V(bp.methyl.mirna)$name,"\\|",2)
V(bp.methyl.mirna)$name = extract_string(V(bp.methyl.mirna)$name,"\\|",1)

vertex_attr(bp.methyl.mirna, "native_nom") = ifelse(V(bp.methyl.mirna)$node_type == "mirna", "miRID", "Symbol")
id = which(V(bp.methyl.mirna)$node_type == "mirna")
vertex_attr(bp.methyl.mirna, "miRID", id) = V(bp.methyl.mirna)$name[id]
id = which(V(bp.methyl.mirna)$node_type == "methyl")
vertex_attr(bp.methyl.mirna, "Symbol", id) = V(bp.methyl.mirna)$name[id]
id = which(V(bp.methyl.mirna)$node_type == "methyl" & V(bp.methyl.mirna)$name %in% rownames(methyl.entrez.dat))
vertex_attr(bp.methyl.mirna, 'EntrezID', id) = methyl.entrez.dat$EntrezID[match(V(bp.methyl.mirna)$name[id], rownames(methyl.entrez.dat))]

#==================================#
# Create data list object for AMEND
#     from DE Analysis Results
#==================================#
# Names of data vectors must correspond to 'name' v.attr of graph
# Add column of names that correspond to node names in graphs
de.res$mrna$node.names = de.res$mrna$EntrezID
de.res$mirna$node.names = rownames(de.res$mirna)
de.res$methyl$node.names = rownames(de.res$methyl)

param.id = list(c(1,1), c(1,2), c(1,3), c(2,1), c(3,1))[[1]] 
data.type = c('logFC', 'P.Value', 'adj.P.Val')[param.id[1]]
pval.weight = c(NULL, 'P.Value', 'adj.P.Val')[param.id[2]]

get_seeds = function(data, data.type, pval.weight=NULL, node.id.col){
  id = which(!is.na(data[, node.id.col]))
  if(data.type == 'logFC'){
    if(!is.null(pval.weight)){
      tmp = data[id, data.type] * (1 - data[id, pval.weight])
    }else tmp = data[id, data.type]
  }
  if(data.type %in% c('P.Value', 'adj.P.Val')){
    tmp = data[,data.type]
  }
  names(tmp) = data[id, node.id.col]
  tmp
} 

data.list = lapply(de.res, get_seeds, data.type = data.type, pval.weight = pval.weight, node.id.col = 'node.names')

#===================================#
# Create BRW list object for AMEND
#     from Survival Analysis Results
#===================================#
# Names of vectors must correspond to 'name' v.attr of graph
# Add column of names that correspond to names in 'hr' (hazard ratios object)
de.res = lapply(de.res, function(x) {
  x$hr = rownames(x)
  x
  })

get_brw_attr = function(brw.data, data.type, pval.weight = FALSE, de.data, node.id.col, brw.col){
  id = which(!is.na(de.data[, node.id.col]))
  tmp.dat = brw.data[de.data[id, brw.col]]
  if(data.type == 'hr'){
    tmp = unlist(lapply(tmp.dat, function(x) x['hr']))
  }
  if(data.type == '1/hr'){
    tmp = unlist(lapply(tmp.dat, function(x) 1/x['hr']))
  }
  if(data.type == '|log(hr)|'){
    tmp = unlist(lapply(tmp.dat, function(x) abs(log(x['hr'])) ))
  }
  if(pval.weight){
    p = unlist(lapply(tmp.dat, function(x) x['pval'] ))
    tmp = tmp * (1 - p)
  }
  names(tmp) = de.data[id, node.id.col]
  tmp
} 

data.type = c('hr', '1/hr', '|log(hr)|')[1]
pval.weight = c(TRUE, FALSE)[2]

brw.list = mapply(FUN = get_brw_attr, brw.data = hr, data.type = data.type, pval.weight = pval.weight, de.data = de.res, node.id.col = 'node.names', brw.col = 'hr')

# lapply(brw.list, summary)

#======#
# AMEND ----
#======#
# AMEND functions
source(paste0(path.amend, "run_AMEND.R"))
source(paste0(path.amend, "create_integrated_graph.R"))
source(paste0(path.amend, "utils.R"))
source(paste0(path.amend, "RandomWalk.R"))
source(paste0(path.amend, "heinz.R"))

input_graphs = list('mirna;mrna' = bp.mirna.mrna,
                   'mirna;methyl' = bp.methyl.mirna,
                   'methyl;mrna' = bp.methyl.mrna,
                   'mrna_string' = g.string,
                   'mrna_path' = g.path,
                   'mrna_aracne' = g.aracne,
                   'mrna_ohmnet' = g.ohmnet)

# Jump & Switch parameters
jp = c(mirna = 1,
       mrna = 0.5,
       methyl = 1)
sp = list(mrna = c(string = 0.5, 
                   path = 0.5, 
                   aracne = 0.5, 
                   ohmnet = 0.5))

# Component and layer seed weights
nw = c(mirna = 0.4, 
       mrna = 0.15,
       methyl = 0.45) 
lw = list(mrna = c(string = 0.25,
                   path = 0.15,
                   aracne = 0.3,
                   ohmnet = 0.3))

id=1
seed.function = c('exp', 'p_value')[id]
fun.params = list(list(DOI = 0), NULL)[[id]]
N = 150
agg.layer = c('string', 'aracne', 'ohmnet', 'path')[1]
aggregate.multiplex = list(primary = paste0("mrna_", agg.layer), agg.method = "mean")


# Create integrated and aggregated graphs
full.graph = create_integrated_graph(graph = input_graphs, data = data.list, brw.attr = brw.list, 
                                     multiplex = TRUE, heterogeneous = TRUE, FUN = seed.function, FUN.params = fun.params, lcc = TRUE)

agg.graph = create_aggregated_graph(graph = full.graph, control = aggregate.multiplex)


#=== Running AMEND ===#

db.adj.methods = list(NULL,
                      list(component = paste0('mrna_', agg.layer), method = 'BS'),
                      list(component = paste0('mrna_', agg.layer), method = 'IN'),
                      list(component = paste0('mrna_', agg.layer), method = 'SDS'))

n_cores = 16

# Non-aggregated
if(1){
  # BRW, non-aggregated
  for(i in seq_along(db.adj.methods)){
    subnet = run_AMEND(graph = input_graphs, n = N, data = data.list, brw.attr = brw.list, aggregate.multiplex = NULL, degree.bias = db.adj.methods[[i]],
                       multiplex = TRUE, heterogeneous = TRUE, normalize = "degree", FUN = seed.function, FUN.params = fun.params,
                       jump.prob = jp, net.weight = nw, switch.layer.prob = sp, layer.weight = lw, verbose = TRUE, in.parallel = TRUE, n.cores = n_cores)
    
    saveRDS(subnet, file = paste0(path.amend.res, "TCGA_n", N, "_db.", ifelse(is.null(db.adj.methods[[i]]), 'none', db.adj.methods[[i]]$method), "_brw.rds"))
  }
  
  # No BRW, non-aggregated
  for(i in seq_along(db.adj.methods)){
    subnet = run_AMEND(graph = input_graphs, n = N, data = data.list, brw.attr = NULL, aggregate.multiplex = NULL, degree.bias = db.adj.methods[[i]],
                       multiplex = TRUE, heterogeneous = TRUE, normalize = "degree", FUN = seed.function, FUN.params = fun.params,
                       jump.prob = jp, net.weight = nw, switch.layer.prob = sp, layer.weight = lw, verbose = TRUE, in.parallel = TRUE, n.cores = n_cores)
    
    saveRDS(subnet, file = paste0(path.amend.res, "TCGA_n", N, "_db.", ifelse(is.null(db.adj.methods[[i]]), 'none', db.adj.methods[[i]]$method), ".rds"))
  }
}

# Aggregated multiplex
if(1){
  # BRW, aggregated multiplex
  for(i in seq_along(db.adj.methods)){
    subnet = run_AMEND(graph = input_graphs, n = N, data = data.list, brw.attr = brw.list, aggregate.multiplex = aggregate.multiplex, degree.bias = db.adj.methods[[i]],
                       multiplex = TRUE, heterogeneous = TRUE, normalize = "degree", FUN = seed.function, FUN.params = fun.params,
                       jump.prob = jp, net.weight = nw, switch.layer.prob = sp, layer.weight = lw, verbose = TRUE, in.parallel = TRUE, n.cores = n_cores)
    
    saveRDS(subnet, file = paste0(path.amend.res, "TCGA_n", N, "_db.", ifelse(is.null(db.adj.methods[[i]]), 'none', db.adj.methods[[i]]$method), "_agg.", agg.layer, "_brw.rds"))
  }
  # No BRW, aggregated multiplex
  for(i in seq_along(db.adj.methods)){
    subnet = run_AMEND(graph = input_graphs, n = N, data = data.list, brw.attr = NULL, aggregate.multiplex = aggregate.multiplex, degree.bias = db.adj.methods[[i]],
                       multiplex = TRUE, heterogeneous = TRUE, normalize = "degree", FUN = seed.function, FUN.params = fun.params,
                       jump.prob = jp, net.weight = nw, switch.layer.prob = sp, layer.weight = lw, verbose = TRUE, in.parallel = TRUE, n.cores = n_cores)
    
    saveRDS(subnet, file = paste0(path.amend.res, "TCGA_n", N, "_db.", ifelse(is.null(db.adj.methods[[i]]), 'none', db.adj.methods[[i]]$method), "_agg.", agg.layer, ".rds"))
  }
}


# '.rds' : nw(mirna)=0.4, nw(methyl)=0.4, nw(mrna)=0.2 <--- Scenario 1
# '_2.rds' : nw(mirna)=0.4, nw(methyl)=0.45, nw(mrna)=0.15 <--- Scenario 2

#=====================#
# Interpreting Results ----
#=====================#
if(0){
  path.end = c("TCGA_n150_db.IN_agg.string_brw.rds", "TCGA_n150_db.IN_agg.string.rds", 'TCGA_n150_db.none_agg.string_brw.rds')
  subnets = vector('list', length(path.end))
  names(subnets) = c('db.brw.agg', 'db.agg', 'brw.agg')
  for(i in seq_along(subnets)){
    file.path = paste0(path.amend.res, path.end[i])
    subnets[[i]] = readRDS(file.path)
    if(!grepl('brw', names(subnets)[i])){
      # Add BRW values to these subnets
      comps = extract_string(extract_string(V(subnets[[i]]$module)$name, "\\|", 2), "_", 1)
      nm = extract_string(V(subnets[[i]]$module)$name, "\\|", 1)
      brw.tmp = numeric(length(nm))
      for(j in seq_along(brw.tmp)) brw.tmp[j] = brw.list[[comps[j]]][match(nm[j], names(brw.list[[comps[j]]]))]
      V(subnets[[i]]$module)$brw.values = brw.tmp
    }
  }

  lapply(subnets, function(x) table(V(x$module)$node_type))
  lapply(subnets, function(x) vcount(x$module))
  lapply(subnets, function(x) ecount(x$module))
  
  # Evaluate effectiveness of BRW. Look at: nodes only in BRW modules... average HR for each module... 
  if(0){
    brw.only.nodes = setdiff(V(subnets$db.brw.agg$module)$name, V(subnets$db.agg$module)$name)
    agg.only.nodes = setdiff(V(subnets$db.agg$module)$name, V(subnets$db.brw.agg$module)$name)
    
    set1 = V(subnets$db.brw.agg$module)$brw.values[V(subnets$db.brw.agg$module)$name %in% brw.only.nodes]
    set2 = V(subnets$db.agg$module)$brw.values[V(subnets$db.agg$module)$name %in% agg.only.nodes]
    
    summary(set1)
    summary(set2)
    
    kruskal.test(x = list(set1, set2))
    
    summary(V(subnets$db.brw.agg$module)$brw.values)
    summary(V(subnets$db.agg$module)$brw.values)
  }
  
  lapply(subnets, function(x) mean(V(x$module)$brw.values))
  
  # Look at average |logFC| for each module
  lapply(subnets,  function(x) mean(abs(V(x$module)$pre.seed.values)))
  
  ## ORA:
  # GO terms as EntrezIDs
  go.desc = readRDS(paste0(path.data, "GO_annotations.rds"))
  go = readRDS(paste0(path.data, "GO_entrez_gene_sets.rds"))
  go.all = do.call(c, go)
  # Merging redundant pathways with "|" character
  if(1){ 
    tmp = unlist(lapply(go.all, function(x) paste(sort(x), collapse = "|")))
    tmp2 = table(tmp)
    tmp3 = vector("list", length(tmp2))
    for(i in seq_along(tmp3)){
      tmp3[[i]] = go.all[[match(names(tmp2)[i], tmp)]]
      nm = names(tmp)[tmp %in% names(tmp2)[i]]
      names(tmp3)[i] = nm[which.min(nchar(nm))]
    }
    go.all = tmp3
  }
  
  #=== Methylation Nodes during AMEND ===#
  methyl.nm = extract_string(V(full.graph)$name[V(full.graph)$node_type == 'methyl'], '\\|', 1)
  res = vector('list', length(subnets))
  for(i in seq_along(subnets)){
    iters = lapply(subnets[[i]]$subnetworks, extract_string, k='\\|', pos=1)
    res[[i]] = numeric(length(iters))
    for(j in seq_along(iters)){
      res[[i]][j] = round(mean(methyl.nm %in% iters[[j]]),3)
    }
  }
  res
  
  
  #=========================#
  # Map EntrezIDs to Symbols ----
  #=========================#
  # uniq.entrez = unique(V(agg.graph)$EntrezID)
  uniq.entrez = unique(unlist(lapply(subnets, function(x) V(x$module)$EntrezID)))
  uniq.entrez = uniq.entrez[!is.na(uniq.entrez)]
  symbol2entrez <- getBM(attributes = c('hgnc_symbol', "entrezgene_id"),
                         filters = 'entrezgene_id',
                         values = uniq.entrez,
                         mart = mart_data)
  # Remove features that didn't map
  symbol2entrez = symbol2entrez[apply(symbol2entrez,1,function(x) all(x != "")), ] 
  subnets = lapply(subnets, function(x){
    id = match(V(x$module)$EntrezID, symbol2entrez[,2])
    V(x$module)$Symbol[!is.na(id)] = symbol2entrez[id[!is.na(id)],1]
    x
  })
  id = match(V(agg.graph)$EntrezID, symbol2entrez[,2])
  V(agg.graph)$Symbol[!is.na(id)] = symbol2entrez[id[!is.na(id)],1]
  
  lapply(subnets, function(x){
    mean(!is.na(V(x$module)$Symbol))
  })
  
  #=== Create 'label' vertex attribute ===#
  subnets = lapply(subnets, function(x){
    V(x$module)$label = ifelse(is.na(V(x$module)$Symbol) & extract_string(V(x$module)$node_type, '_', 1) == 'mrna', V(x$module)$EntrezID, 
                               ifelse(is.na(V(x$module)$Symbol), V(x$module)$miRID, V(x$module)$Symbol))
    V(x$module)$label.extra = ifelse(grepl('_', V(x$module)$node_type), paste(V(x$module)$label, extract_string(V(x$module)$node_type,'_',2), sep='|'), V(x$module)$label)
    x
  })
  lapply(subnets, function(x){
    head(V(x$module)$label)
  })
  
  #=== Looking at non-aggregated subnets ===#
  ids=which(!grepl('agg', names(subnets)))
  res = vector('list', length(ids))
  for(i in seq_along(ids)){
    tmp = table(V(subnets[[ids[i]]]$module)$label)
    res[[i]] = sort(tmp, decreasing = TRUE)
  }
  res
  
  #==================#
  # ORA with GO Terms ----
  #==================#
  library(fgsea)
  
  all.entrezIDs = extract_string(V(full.graph)$name, '\\|', 1)
  all.entrezIDs = unique(all.entrezIDs[grepl('mrna', V(full.graph)$node_type)])
  
  ora = vector('list', length(subnets))
  names(ora) = names(subnets)
  for(i in seq_along(ora)){
    tmp = extract_string(V(subnets[[i]]$module)$name, '\\|', 1)
    tmp = unique(tmp[grepl('mrna', extract_string(V(subnets[[i]]$module)$name, '\\|', 2))])
    # ORA on GO aspects separately
    if(1){
      ora.tmp = NULL
      for(j in seq_along(go)){
        ora.tmp2 = fora(pathways = go[[j]], 
                        genes = tmp, 
                        universe = all.entrezIDs) %>%
          dplyr::filter(padj <= 0.01) %>%
          dplyr::mutate(overlapGenes = unlist(lapply(overlapGenes, function(x) paste(x, collapse = ", "))),
                        GO_ID = pathway,
                        GO_aspect = extract_string(names(go)[j], '_', 2))
        ora.tmp = rbind(ora.tmp, ora.tmp2, fill=TRUE)
      }
      ora.tmp = ora.tmp[order(ora.tmp$padj)]
      ora[[i]] = ora.tmp
    }
    # ORA on all GO aspects simultaneously
    if(0){
      ora[[i]] = fora(pathways = go.all,
                      genes = tmp,
                      universe = all.entrezIDs) %>%
        dplyr::filter(padj <= 0.01) %>%
        dplyr::mutate(overlapGenes = unlist(lapply(overlapGenes, function(x) paste(x, collapse = ", "))))
      
      ora[[i]]$GO_ID = extract_string(ora[[i]]$pathway, '\\.', 2)
      ora[[i]]$GO_aspect = ifelse(extract_string(ora[[i]]$pathway, '\\.', 1) == 'go_p', 'p',
                                  ifelse(extract_string(ora[[i]]$pathway, '\\.', 1) == 'go_f', 'f', 'c'))
    }

    tmp = character(nrow(ora[[i]]))
    for(j in seq_along(tmp)){
      tmp[j] = go.desc[[ora[[i]]$GO_aspect[j]]]$term[rownames(go.desc[[ora[[i]]$GO_aspect[j]]]) == ora[[i]]$GO_ID[j]]
    }
    ora[[i]]$pathway = tmp
  }
  lapply(ora, nrow)
  lapply(ora, head)
  lapply(ora, function(x) table(x$GO_aspect))
  lapply(ora, function(x) x$pathway[x$GO_aspect == 'p'])
  
  #===================#
  # Disease Enrichment ----
  #===================#
  library(DOSE)
  p.adj.cutoff = 0.01
  dgn = do = vector('list', length(subnets))
  names(dgn) = names(do) = names(subnets)
  for(i in seq_along(subnets)){
    message(i)
    tmp = V(subnets[[i]]$module)$EntrezID
    tmp = tmp[!is.na(tmp)]
    all.tmp = unique(V(full.graph)$EntrezID)
    all.tmp = all.tmp[!is.na(all.tmp)]
    dgn[[i]] = enrichDGN(gene = tmp,
                    pAdjustMethod = "BH",
                    universe = all.tmp)
    dgn[[i]]@result = dgn[[i]]@result[dgn[[i]]@result$p.adjust <= p.adj.cutoff,]
    do[[i]] = enrichDO(gene = tmp,
                  ont = "DO",
                  pAdjustMethod = "BH",
                  universe = all.tmp)
    do[[i]]@result = do[[i]]@result[do[[i]]@result$p.adjust <= p.adj.cutoff,]
  }
  
  #=======================#
  # Visualize in Cytoscape ----
  #=======================#
  library(RCy3)
  library(RColorBrewer)
  
  # Style 1
  if(0){
    for(i in seq_along(subnets)){
      message(i)
      createNetworkFromIgraph(subnets[[i]]$module, title = paste0("TCGA_n", subnets[[i]]$input_params$n, '_', names(subnets)[i]), collection = "TCGA")
      # net.id = getNetworkSuid()
      
      ### Set the visual style
      # Looking at color schemes
      if(0){ 
        display.brewer.all()
        display.brewer.pal(n = 9, name = "RdBu")
        brewer.pal(n = 11, name = "RdBu")
      }
      color.theme = "RdBu"
      # The above may need to change depending on if you want divergent or gradient color themes, which will depend on the data type used in AMEND
      
      style.name = names(subnets)[i]
      createVisualStyle(style.name)
      
      # Set node shape and size
      lockNodeDimensions(new.state = TRUE, style.name = style.name)
      setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SHAPE", value = "ellipse"), style.name = style.name)
      setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SIZE", value = 40), style.name = style.name)
      
      # Set node color as a function of Experimental data used in AMEND. This will need to be modified for different data types
      color.scheme = brewer.pal(n = 3, name = color.theme)
      v.attr.name = "pre.seed.values" # The name of the vertex attribute used as seed values in AMEND
      v.range = c(min(vertex_attr(subnets[[i]]$module, v.attr.name)), max(vertex_attr(subnets[[i]]$module, v.attr.name)))
      setNodeColorMapping(table.column = v.attr.name, table.column.values = c(v.range[1], 0, v.range[2]), colors = color.scheme, style.name = style.name)
      # Node labels, position, and font size
      tmp = ifelse(grepl('agg', names(subnets)[i]), 'label', 'label.extra')
      setNodeLabelMapping(table.column = tmp, style.name = style.name)
      setNodeLabelPositionDefault(new.nodeAnchor = "N", new.graphicAnchor = "S", new.justification = "c", new.xOffset = 0, new.yOffset = 0, style.name = style.name)
      setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_SIZE", value = 28), style.name = style.name)
      setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_FACE", value = "AlBayan-Bold"), style.name = style.name)
      # Set edge style and width
      setEdgeLineStyleDefault(new.line.style = "DOT", style.name = style.name)
      setEdgeLineWidthDefault(new.width = 3, style.name = style.name)
      
      setVisualStyle(style.name)
    }
  }
  # Style 2: node.color=seed value, border.color=component, label=name|layer or name|component
  if(0){
    for(i in seq_along(subnets)){
      message(i)
      createNetworkFromIgraph(subnets[[i]]$module, title = paste0("TCGA_n", subnets[[i]]$input_params$n, '_', names(subnets)[i]), collection = "TCGA")
      layoutNetwork(paste("force-directed", "defaultSpringLength=70 defaultSpringCoefficient=0.000007", sep = " "))
      
      ### Set the visual style
      # Looking at color schemes
      if(0){ 
        display.brewer.all()
        display.brewer.pal(n = 9, name = "RdBu")
        brewer.pal(n = 11, name = "RdBu")
      }
      color.theme = "RdBu"
      # The above may need to change depending on if you want divergent or gradient color themes, which will depend on the data type used in AMEND
      
      style.name = paste0(names(subnets)[i])
      createVisualStyle(style.name)
      
      # Set node shape and size
      lockNodeDimensions(new.state = TRUE, style.name = style.name)
      setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SHAPE", value = "ellipse"), style.name = style.name)
      setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SIZE", value = 40), style.name = style.name)
      
      # Set node color as a function of Experimental data used in AMEND. This will need to be modified for different data types
      color.scheme = brewer.pal(n = 3, name = "Greys")
      v.attr.name = "seeds" # The name of the vertex attribute used as seed values in AMEND
      v.range = c(min(vertex_attr(subnets[[i]]$module, v.attr.name)), 
                  median(vertex_attr(subnets[[i]]$module, v.attr.name)), 
                  max(vertex_attr(subnets[[i]]$module, v.attr.name)))
      # v.range = c(0, 0.5, 1)
      setNodeColorMapping(table.column = v.attr.name, table.column.values = v.range, colors = color.scheme, style.name = style.name)
      # Node labels, position, and font size
      tmp = ifelse(grepl('agg', names(subnets)[i]), 'label', 'label.extra2')
      setNodeLabelMapping(table.column = tmp, style.name = style.name)
      setNodeLabelPositionDefault(new.nodeAnchor = "N", new.graphicAnchor = "S", new.justification = "c", new.xOffset = 0, new.yOffset = 0, style.name = style.name)
      setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_SIZE", value = 28), style.name = style.name)
      setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_FACE", value = "AlBayan-Bold"), style.name = style.name)
      # Set edge style and width
      setEdgeLineStyleDefault(new.line.style = "DOT", style.name = style.name)
      setEdgeLineWidthDefault(new.width = 3, style.name = style.name)
      
      # Set node border
      setNodeBorderWidthBypass(node.names = V(subnets[[i]]$module)$name, new.sizes = 9)
      # Set node border color
      color.scheme = brewer.pal(n = 9, name = "Set1")[c(1,2,5,8)] # red, blue, orange, pink
      # !!!!!!!!!!!!!!!!!!!!
      if(grepl('agg', names(subnets)[i])){
        nt = extract_string(V(subnets[[i]]$module)$node_type, '_', 1)
        new.node.colors = ifelse(nt == "mrna", color.scheme[1], ifelse(nt == "methyl", color.scheme[3], color.scheme[2]))
        # nt = extract_string(V(subnets[[i]]$module)$name, "\\|", 2)
        # new.node.colors = ifelse(nt != "meta", color.scheme[1], color.scheme[2])
      }else{
        nt = extract_string(V(subnets[[i]]$module)$node_type, '_', 1)
        new.node.colors = ifelse(nt == "t", color.scheme[1], ifelse(nt == "p", color.scheme[3], ifelse(nt == "ph", color.scheme[4], color.scheme[2])))
      }
      setNodeBorderColorBypass(node.names = V(subnets[[i]]$module)$name, new.colors = new.node.colors)
      
      setVisualStyle(style.name)
    }
  }
  
  
  #====================================#
  # Visualize ORA Pathways in Cytoscape ----
  #====================================#
  library(RCy3)
  library(RColorBrewer)
  
  subnet = subnets$db.brw.agg
  ora.paths = ora[[1]]
  ora.paths = ora.paths[order(ora.paths$padj),]
  pathway.list = go.all
  # Change names of pathway.list
  if(1){
    go.aspect = extract_string(extract_string(names(go.all), '\\.',1), '_', 2)
    go.id = extract_string(names(go.all), '\\.',2)
    tmp = character(length(go.all))
    for(j in seq_along(tmp)){
      tmp[j] = go.desc[[go.aspect[j]]]$term[rownames(go.desc[[go.aspect[j]]]) == go.id[j]]
    }
    names(pathway.list) = tmp
  }

  nested_pathways = function(enriched.pathways, rep.path.ids, all.pathways, threshold = 0.8, first.n = NULL, parallel = FALSE, nclust, padj = TRUE){
    require(doParallel)
    require(foreach)
    
    if(nrow(enriched.pathways) == 0){
      stop("nrow of \'enriched.pathways\' is zero")
    }else if(nrow(enriched.pathways) == 1){
      ora_path = enriched.pathways$pathway
      g = igraph::make_empty_graph(n = length(ora_path))
      # vlabs = str_replace(ora_path, "Mus musculus: ", "") 
      vlabs = ora_path
      V(g)$name = vlabs
      temp = numeric(vcount(g))
      for(i in seq_along(temp)){
        temp[i] = length(all.pathways[[ora_path[i]]])
      }
      V(g)$size = temp
      V(g)$pval = ifelse(padj, enriched.pathways$padj, enriched.pathways$pval)
      col.tmp = rep(0, igraph::vcount(g))
      col.tmp[rep.path.ids] = seq_along(rep.path.ids)
      V(g)$color = col.tmp
      return(g)
    }
    
    if(!is.null(first.n)) enriched.pathways = enriched.pathways[1:min(first.n, nrow(enriched.pathways)),]
    
    nested_index = function(a, b){
      intersection = length(intersect(a, b))
      denom = min(length(a), length(b))
      return(intersection / denom)
    }
    
    ora_path = enriched.pathways$pathway
    p.id1 = c()
    p.id2 = c()
    sim.mat = matrix(0, nrow = length(ora_path), ncol = length(ora_path), dimnames = list(ora_path, ora_path))
    for(i in 1:(length(ora_path)-1)){
      for(j in (i+1):length(ora_path)){
        sim.mat[i,j] = nested_index(all.pathways[[ora_path[i]]], all.pathways[[ora_path[j]]])
        sim.mat[j,i] = sim.mat[i,j]
        if(sim.mat[i,j] > threshold){ 
          p.id1 = c(p.id1, ifelse(min(length(all.pathways[[ora_path[i]]]), length(all.pathways[[ora_path[j]]])) == length(all.pathways[[ora_path[i]]]), i, j))
          p.id2 = c(p.id2, ifelse(min(length(all.pathways[[ora_path[i]]]), length(all.pathways[[ora_path[j]]])) == length(all.pathways[[ora_path[i]]]), j, i))
        }
      }
    }
    
    # Finding nested pathways 
    # nest.mat is an edgelist, where a vertex in col 1 is nested within a vertex of col 2 (vertex = pathway)
    if(is.null(p.id1) || is.null(p.id2)){
      g = igraph::make_empty_graph(n = length(ora_path))
      # vlabs = str_replace(ora_path, "Mus musculus: ", "") 
      vlabs = ora_path
      V(g)$name = vlabs
      temp = numeric(vcount(g))
      for(i in seq_along(temp)){
        temp[i] = length(all.pathways[[ora_path[i]]])
      }
      V(g)$size = temp
      V(g)$pval = ifelse(padj, enriched.pathways$padj, enriched.pathways$pval)
      col.tmp = rep(0, igraph::vcount(g))
      col.tmp[rep.path.ids] = seq_along(rep.path.ids)
      V(g)$color = col.tmp
      return(g)
    }else{
      nest.mat = matrix(c(p.id1, p.id2), ncol = 2)
      g = igraph::graph_from_edgelist(nest.mat, directed = TRUE)
      for(i in seq_along(setdiff(1:length(ora_path), V(g)))) g = igraph::add_vertices(g, 1)
    }
    
    # Getting longest simple, directed paths between pair of vertices in edgelist
    # Rationale: If a pathway is nested w/in a larger pathway, it may be nested in other pathways that are themselves nested in that larger pathway.
    #   I want to see the hierarchy of nested-ness
    if(parallel){ # parallel FOR loop
      cl = makeForkCluster(nclust, outfile = "")
      registerDoParallel(cl)
      nest.id0 = foreach(i = 1:nrow(nest.mat), .packages = "igraph") %dopar% {
        message(paste0(i, ":", nrow(nest.mat)))
        asp = all_simple_paths(graph = g, from = nest.mat[i,1], to = nest.mat[i,2], mode = "out")
        asp[[which.max(do.call("c", lapply(asp, length)))]]
      }
      stopCluster(cl)
    }else{ # normal FOR loop
      nest.id0 = vector("list", nrow(nest.mat))
      for(i in 1:nrow(nest.mat)){ # This for loop may take a long time
        message(paste0(i, ":", nrow(nest.mat)))
        asp = all_simple_paths(graph = g, from = nest.mat[i,1], to = nest.mat[i,2], mode = "out")
        nest.id0[[i]] = asp[[which.max(do.call("c", lapply(asp, length)))]]
        # if(i %% 50 == 0) message(paste0(i, ":", nrow(nest.mat)))
      }
    }
    
    # This gets rid of redundant paths
    nest.id = list()
    for(i in 1:length(nest.id0)){
      cond = T
      for(j in (1:length(nest.id0))[(1:length(nest.id0)) != i]){
        if(all(nest.id0[[i]] %in% nest.id0[[j]])) cond = F
      }
      if(cond) nest.id[[length(nest.id) + 1]] = nest.id0[[i]]
    }
    
    # Converting vertex ids to pathway names
    nested.paths = vector(mode = "list", length = length(nest.id))
    for(i in 1:length(nest.id)){
      nested.paths[[i]] = ora_path[nest.id[[i]]]
    }
    # Going from most nested to least nested (smallest to largest)
    
    # Getting graph showing nested hierarchy structure
    el = matrix(nrow = 1000, ncol = 2)
    r.id = 1
    for(i in seq_along(nest.id)){
      x = as.numeric(nest.id[[i]])
      for(j in 1:(length(x) - 1)){
        new.row = x[j:(j+1)]
        if(any(apply(el, 1, function(y) all(new.row %in% y)))) next
        el[r.id,] = new.row 
        r.id = r.id + 1
      }
    }
    el = na.omit(el)
    g = igraph::graph_from_edgelist(el = el, directed = TRUE)
    for(i in seq_along(setdiff(1:length(ora_path), V(g)))) g = igraph::add_vertices(g, 1)
    # vlabs = str_replace(ora_path, "Mus musculus: ", "") 
    vlabs = ora_path
    V(g)$name = vlabs
    temp = numeric(vcount(g))
    for(i in seq_along(temp)){
      temp[i] = length(all.pathways[[ora_path[i]]])
    }
    V(g)$size = temp
    V(g)$pval = enriched.pathways$padj
    col.tmp = rep(0, igraph::vcount(g))
    col.tmp[rep.path.ids] = seq_along(rep.path.ids)
    V(g)$color = col.tmp
    return(g)
  }
  
  # Get that pathways that each node in module is associated with
  node.pathway = vector("list", vcount(subnet$module)); names(node.pathway) = V(subnet$module)$node_label
  for(i in seq_along(node.pathway)){
    node.pathway[[i]] = ora.paths$pathway[unlist(lapply(strsplit(ora.paths$overlapGenes, ", "), function(x) V(subnet$module)$EntrezID[i] %in% x))]
  }
  
  ## Getting cellular functions associated with each cluster in module
  set.seed(2088) # Louvain is stochastic, so we set a seed value for reproducibility 
  # Running Louvain clustering algorithm on AMEND module to get clusters
  cl = igraph::cluster_louvain(graph = subnet$module, resolution = 1)
  cl.id = unique(cl$membership)
  module.cluster = character(length(cl.id))
  # For each cluster, identify pathway that contains the most nodes of that cluster
  for(j in seq_along(cl.id)){
    node.names = cl$names[which(cl$membership == cl.id[j])]
    tmp = unlist(node.pathway[which(V(subnet$module)$name %in% node.names)])
    if(length(tmp) == 0){ # If the selected nodes aren't part of any pathway...
      module.cluster[j] = "unknown"
    }else{
      p.count = table(tmp)
      tmp.id = which(p.count == max(p.count))
      # In case of ties for most occurring pathway, take largest pathway
      if(length(tmp.id) > 1) tmp.id = tmp.id[which.max(unlist(lapply(pathway.list[names(p.count)[tmp.id]], length)))]
      module.cluster[j] = names(p.count)[tmp.id]
    }
    names(module.cluster)[j] = paste(V(subnet$module)$label[V(subnet$module)$name %in% node.names], collapse = ", ")
  }
  if(any(table(module.cluster) > 1)){ # Merge clusters that have same Representative Pathway
    tmp = table(module.cluster)
    tmp.id = which(tmp > 1)
    for(j in seq_along(tmp.id)){
      tmp.id2 = which(module.cluster %in% names(tmp)[tmp.id[j]])
      keep.id = tmp.id2[1]
      rm.id = tmp.id2[-1]
      names(module.cluster)[keep.id] = paste(c(names(module.cluster)[keep.id], names(module.cluster)[rm.id]), collapse = ", ")
      module.cluster = module.cluster[-rm.id]
    }
  }
  tmp = data.frame(.function = module.cluster, genes = names(module.cluster))
  tmp = tmp[tmp$.function != "unknown",]
  module.cluster = tmp
  
  ## Visualize pathways in Cytoscape
  top.k.pathways = 30
  threshold = 0.6
  tkp = min(top.k.pathways, nrow(ora.paths))
  p.tmp = ora.paths[1:tkp,]
  rp.tmp = module.cluster$.function[!module.cluster$.function %in% p.tmp$pathway]
  if(length(rp.tmp) == 0){
    rp.ids = which(p.tmp$pathway %in% module.cluster$.function)
  }else{
    id.tmp = which(!p.tmp$pathway %in% module.cluster$.function)
    p.tmp = rbind2(p.tmp[-tail(id.tmp, length(rp.tmp)),], ora.paths[ora.paths$pathway %in% rp.tmp,])
    rp.ids = which(p.tmp$pathway %in% module.cluster$.function)
  }
  g.path = nested_pathways(enriched.pathways = p.tmp, rep.path.ids = rp.ids, all.pathways = pathway.list, threshold = threshold, first.n = NULL, parallel = TRUE, nclust = 4)
  createNetworkFromIgraph(g.path, title = "TCGA DAG", collection = "TCGA DAG")
  layoutNetwork("hierarchical")
  
  ## Set Visual Style
  style.name = "tcga_dag"
  createVisualStyle(style.name)
  
  # Set node shape and size
  lockNodeDimensions(new.state = TRUE, style.name = style.name)
  setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SHAPE", value = "ellipse"), style.name = style.name)
  setNodeSizeMapping(table.column = "size", table.column.values = c(10, 500, 1000), 
                     sizes = c(15, 55, 100), style.name = style.name)
  
  # Set node color
  gray = brewer.pal(n = 9, name = "Set1")[9]
  col.tmp = rep(gray, igraph::vcount(g.path))
  if(sum(V(g.path)$color != 0) > 8 & sum(V(g.path)$color != 0) <= 16){
    col.tmp[V(g.path)$color != 0][1:8] = brewer.pal(n = 8, name = "Set1")[1:8]
    col.tmp[V(g.path)$color != 0][9:sum(V(g.path)$color != 0)] = brewer.pal(n = 8, name = "Set2")[1:(sum(V(g.path)$color != 0) - 8)]
  }else if(sum(V(g.path)$color != 0) <= 8){
    col.tmp[V(g.path)$color != 0] = brewer.pal(n = max(sum(V(g.path)$color != 0),3), name = "Set1")[1:sum(V(g.path)$color != 0)]
  }else stop('too many clusters')
  
  setNodeColorMapping(table.column = "color", table.column.values = V(g.path)$color, colors = col.tmp, mapping.type = "d", style.name = style.name)
  
  # Node labels, position, and font size
  setNodeLabelMapping(table.column = "name", style.name = style.name)
  setNodeLabelPositionDefault(new.nodeAnchor = "N", new.graphicAnchor = "S", new.justification = "c", new.xOffset = 0, new.yOffset = 0, style.name = style.name)
  setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_SIZE", value = 20), style.name = style.name)
  font_style = c("Rockwell-Bold", "ArialNarrow-Bold", "Arial-BoldMT", "Arial-Black")[4]
  setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_FACE", value = font_style), style.name = style.name)
  # Set edge style and width
  setEdgeLineStyleDefault(new.line.style = "DOT", style.name = style.name)
  setEdgeLineWidthDefault(new.width = 4, style.name = style.name)
  setEdgeTargetArrowShapeDefault(new.shape = "DELTA", style.name = style.name)
  setVisualStyle(style.name)
  
  ### Select Nodes in AMEND module corresponding to Representative pathways
  module.cluster$.function # Representative pathways
  n.cl = 8
  module.cluster$.function[n.cl]
  # Select nodes that are in the cluster corresponding to this representative pathway
  selectNodes(nodes = unlist(strsplit(module.cluster$genes[n.cl], ", ")), by.col = "label",
              preserve.current.selection = FALSE)
}


#================#
# DIABLO Analysis ----
#================#
library(mixOmics)

set.seed(2707083)

#=== Process TCGA Data ===#
# Keep only samples that appear in all 3 datasets
com = intersect(colnames(omics.dat$mrna), intersect(colnames(omics.dat$mirna), colnames(omics.dat$methyl)))
diablo_dat = lapply(omics.dat, function(x) {
  t(x[,match(com, colnames(x))])
})


# Get vector of sample class memberships
tt_tmp = extract_string(rownames(diablo_dat$mrna), "-", 4)
Y = ifelse(tt_tmp == "01", "tumor", ifelse(tt_tmp == "11", "normal", "other"))

### DIABLO
## Create Design Matrix... Lower values = better class discrimination, Higher values = more correlated features across datasets
design = matrix(0.5, ncol = length(diablo_dat), nrow = length(diablo_dat), 
                dimnames = list(names(diablo_dat), names(diablo_dat)))
diag(design) = 0

## Choose number of components
dist.name = "mahalanobis.dist"
if(0){
  diablo.tcga <- block.plsda(diablo_dat, Y, ncomp = 5, design = design)
  
  perf.diablo.tcga = perf(diablo.tcga, validation = 'Mfold', folds = 10, nrepeat = 10)
  
  perf.diablo.tcga$error.rate  # Lists the different types of error rates
  
  # Plot of the error rates based on weighted vote
  plot(perf.diablo.tcga)
  
  perf.diablo.tcga$choice.ncomp$WeightedVote
  
  ncomp <- perf.diablo.tcga$choice.ncomp$WeightedVote["Overall.BER", dist.name]
}

ncomp = 5

# Number of variables to select for each dataset
if(0){
  test.keepX <- list(mrna = c(seq(5, 30, 5)),
                     mirna = c(seq(5, 30, 5)),
                     methyl = c(seq(5, 30, 5)))
  tune.diablo.tcga <- tune.block.splsda(diablo_dat, Y, ncomp = ncomp, 
                                        test.keepX = test.keepX, design = design,
                                        validation = 'Mfold', folds = 10, nrepeat = 1, 
                                        BPPARAM = BiocParallel::SnowParam(workers = 4),
                                        dist = dist.name)
  list.keepX <- tune.diablo.tcga$choice.keepX
  # Recommended
  list.keepX = list(mrna = c(5, 5, 5, 5, 5),
                    mirna = c(5, 5, 5, 5, 5),
                    methyl = c(5, 5, 5, 5, 5)) 
}

list.keepX = list(mrna = c(150, 5, 5, 5, 5),
                  mirna = c(6, 5, 5, 5, 5),
                  methyl = c(5, 5, 5, 5, 5)) 


# Fit Final Model
diablo.tcga <- block.splsda(diablo_dat, Y, ncomp = ncomp, 
                            keepX = list.keepX, design = design)

factor_id = 1
mrna.d = selectVar(diablo.tcga, block = 'mrna', comp = factor_id)$mrna$name
mirna.d = selectVar(diablo.tcga, block = 'mirna', comp = factor_id)$mirna$name
methyl.d =selectVar(diablo.tcga, block = 'methyl', comp = factor_id)$methyl$name


# Load in final AMEND module (n=150, using IN degree bias adjustment, hazard ratios as biased random walk attribute, and aggregating mRNA component using STRING edges)
amend = readRDS(paste0(path.amend.res, "TCGA_n150_db.IN_agg.string_brw.rds"))

## Mapping Entrez IDs to Symbols
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://jan2024.archive.ensembl.org")
mart_data <- useDataset("hsapiens_gene_ensembl", mart = mart)
uniq.entrez = unique(V(amend$module)$EntrezID)
uniq.entrez = uniq.entrez[!is.na(uniq.entrez)]
symbol2entrez <- getBM(attributes = c('hgnc_symbol', "entrezgene_id"),
                       filters = 'entrezgene_id',
                       values = uniq.entrez,
                       mart = mart_data)
# Remove features that didn't map
symbol2entrez = symbol2entrez[apply(symbol2entrez,1,function(x) all(x != "")), ] 
id = match(V(amend$module)$EntrezID, symbol2entrez[,2])
V(amend$module)$Symbol[!is.na(id)] = symbol2entrez[id[!is.na(id)],1]


## mRNA
mrna.a = V(amend$module)$Symbol[extract_string(V(amend$module)$name, "\\|", 2) == "mrna"]
mrna.ji = length(intersect(mrna.a, mrna.d)) / min(length(mrna.a), length(mrna.d)) # Jaccard Index: AMEND & DIABLO

## miRNA
mirna.a = extract_string(V(amend$module)$name[extract_string(V(amend$module)$name, "\\|", 2) == "mirna"], "\\|", 1)
mirna.ji = length(intersect(mirna.a, mirna.d)) / min(length(mirna.a), length(mirna.d)) # Jaccard Index: AMEND & DIABLO


# AMEND & DIABLO... Figure 8
### LogFC weighted by 1 minus P-value
par(mfrow = c(1,2))
## mRNA
logfc.a = abs(de.res$mrna$logFC[rownames(de.res$mrna) %in% mrna.a])
logfc.d = abs(de.res$mrna$logFC[rownames(de.res$mrna) %in% mrna.d])
pval.a = de.res$mrna$P.Value[rownames(de.res$mrna) %in% mrna.a]
pval.d = de.res$mrna$P.Value[rownames(de.res$mrna) %in% mrna.d]
score.a = logfc.a * (1 - pval.a)
score.d = logfc.d * (1 - pval.d)
dat = data.frame(score = c(score.a, score.d), method = c(rep("AMEND", length(score.a)), rep("DIABLO", length(score.d))))
# kw_pval = kruskal.test(score ~ method, data = dat)$p.value
w_pval = wilcox.test(score ~ method, data = dat)$p.value
boxplot(score ~ method, data = dat, xlab = "Method", ylab = "|logFC| * (1 - pvalue)", 
        main = paste0("mRNA Features (MW p-value", ifelse(w_pval < 1e-4, " < 0.0001)", paste0(" = ", round(w_pval, 4), ")"))))

## miRNA
logfc.a = abs(de.res$mirna$logFC[rownames(de.res$mirna) %in% mirna.a])
logfc.d = abs(de.res$mirna$logFC[rownames(de.res$mirna) %in% mirna.d])
pval.a = de.res$mirna$P.Value[rownames(de.res$mirna) %in% mirna.a]
pval.d = de.res$mirna$P.Value[rownames(de.res$mirna) %in% mirna.d]
score.a = logfc.a * (1 - pval.a)
score.d = logfc.d * (1 - pval.d)
dat = data.frame(score = c(score.a, score.d), method = c(rep("AMEND", length(score.a)), rep("DIABLO", length(score.d))))
# kw_pval = kruskal.test(score ~ method, data = dat)$p.value
w_pval = wilcox.test(score ~ method, data = dat)$p.value
stripchart(score ~ method, data = dat, vertical = TRUE, method = "overplot",
           pch = 7, col = "blue", xlab = "Method", ylab = "|logFC| * (1 - pvalue)",
           main = paste0("miRNA Features (MW p-value", ifelse(w_pval < 1e-4, " < 0.0001)", paste0(" = ", round(w_pval, 4), ")"))), 
           at = c(1.25, 1.75), cex = 1.5)