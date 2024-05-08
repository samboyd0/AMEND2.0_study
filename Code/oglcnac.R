#=======================#
# O-GlcNAc Study:
#=======================#
# ECI type: OGT-KO 1wk vs. 2wk (+)
# Networks: NISE PPIN, PPIN from STRING, MMIN from STITCH

rm(list = ls())

# Set path names
path.data = "/path/to/data/"
path.amend = "/path/to/AMEND_code/"
path.amend.res = "/path/to/results/"

# Attach libraries
library(igraph)

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

#====================#
# DE Analysis Results ----
#====================#
ECI = readRDS(paste0(path.data, 'ECI.rds'))

# Getting a list of names (with '.' stripped off 
eci.symbols = lapply(ECI[1], function(x) extract_string(names(x), '\\.', 1))
eci.acc = lapply(ECI[2:3], names)

#=====================#
# Interaction Networks ----
#=====================#
# 3 components representing T, P, & Ph, each with 2 layers: STRING physical and NISE
# 1 component representing M
# 4 bipartite graphs: T-P, T-Ph, P-Ph, M-(T/P/Ph)

ppi_edge_threshold = 0.7

#=== STRING PPIN 
ppi = data.table::fread(file = paste0(path.data, "10090.protein.links.v12.0.txt.gz"), header = T) %>%
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

g.ppi = igraph::graph_from_edgelist(el = ppi[,1:2], directed = FALSE)
E(g.ppi)$weight = as.numeric(ppi[,3])
g.ppi = largest_connected_component(g.ppi)

V(g.ppi)$native_nom = 'ENSPID'
V(g.ppi)$ENSPID = V(g.ppi)$name

#========#
# Mapping ----
#========#
# list the databases available from biomaRt and archived versions
if(0){
  listMarts() 
  listEnsemblArchives()
} 
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://jan2024.archive.ensembl.org")
# List datasets available within a certain database
if(0) listDatasets(mart) 
mart_data <- useDataset("mmusculus_gene_ensembl", mart = mart)
# Check the names of available filters & attributes
if(0){
  View(listFilters(mart_data))
  View(listAttributes(mart_data))  
}
# Attributes are what you want returned from your mapping (i.e., what you are mapping to)
# Filters are what you are trying to map 
## Symbols to Ensembl peptide IDs
uniq.symbols = unique(unlist(eci.symbols))
symbol2ensp <- getBM(attributes = c("mgi_symbol", "ensembl_peptide_id"),
                     filters = "mgi_symbol",
                     values = uniq.symbols,
                     mart = mart_data)
# Remove features that didn't map
symbol2ensp = symbol2ensp[apply(symbol2ensp,1,function(x) all(x != "")), ]

uniq.acc = unique(unlist(eci.acc))
acc2ensp = getBM(attributes = c("refseq_peptide", "ensembl_peptide_id"),
                 filters = "refseq_peptide",
                 values = uniq.acc,
                 mart = mart_data)

# Remove features that didn't map
acc2ensp = acc2ensp[apply(acc2ensp,1,function(x) all(x != "")), ]

#=== STRING PPIN Construction ----
# - take largest connected component of ppi
# - Get mapping between ENSPID and symbols/ensg/accession IDs
# - Get KNN graph. Increase knn until connected.

knn = 1
g.strings = vector('list', 3)
names(g.strings) = c('t', 'p', 'ph')
for(i in seq_along(g.strings)){
  if(names(g.strings)[i] == 't'){
    id = which(V(g.ppi)$name %in% symbol2ensp[symbol2ensp[,2] %in% eci.symbols[[names(g.strings)[i]]],1])
  }else{
    id = which(V(g.ppi)$name %in% acc2ensp[acc2ensp[,1] %in% eci.acc[[names(g.strings)[i]]], 2])
  }
  t.ids = unique(as.numeric(unlist(igraph::neighborhood(graph = g.ppi, order = knn, nodes = id))))
  tmp = igraph::induced_subgraph(g.ppi, t.ids)
  if(names(g.strings)[i] == 't'){
    id = match(V(tmp)$name, symbol2ensp[,1])
    vertex_attr(tmp, 'Symbol', which(!is.na(id))) = symbol2ensp[id[!is.na(id)],2]
  }else{
    id = match(V(tmp)$name, acc2ensp[,2])
    vertex_attr(tmp, 'RefSeq', which(!is.na(id))) = acc2ensp[id[!is.na(id)],1]
  }
  V(tmp)$node_type = paste(names(g.strings)[i], 'string', sep='_')
  g.strings[[i]] = largest_connected_component(tmp)
}
# unname(unlist(lapply(g.strings, vcount)))
# 14718 13967 12643 (knn=2)
# 14085 10198  5736 (knn=1)
#==== NISE PPIN Construction ----
# - NISE is much smaller than STRING (only ~820 nodes).
# - Keep all nodes.
nise = readRDS(paste0(path.data, "NISE_graph.rds"))

# Remove unnecessary vertex attributes
v = vertex_attr_names(nise)
v = setdiff(v, c('name', 'native_nom', 'Symbol', 'ENSPID'))
for(i in seq_along(v)) nise = igraph::delete_vertex_attr(nise, v[i])

g.nise = vector('list', length(g.strings))
names(g.nise) = names(g.strings)
for(i in seq_along(g.nise)){
  g.nise[[i]] = nise
  V(g.nise[[i]])$node_type = paste(names(g.nise)[i], 'nise', sep='_')
}

#=== MMI Network of merged isomers ===#
mmi_edge_threshold = 0.9

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

g.meta = graph_from_edgelist(el = mmi[,1:2], directed = FALSE)
E(g.meta)$weight = as.numeric(mmi[,3])
g.meta = largest_connected_component(g.meta)
V(g.meta)$node_type = 'meta'
V(g.meta)$native_nom = 'PubChemID'
V(g.meta)$PubChemID = V(g.meta)$name

# Get KNN graph of nodes that have experimental data
knn = 1
id = which(V(g.meta)$name %in% names(ECI$meta))
t.ids = unique(as.numeric(unlist(igraph::neighborhood(graph = g.meta, order = knn, nodes = id))))
tmp = igraph::induced_subgraph(g.meta, t.ids)
g.meta = largest_connected_component(tmp)

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

bp[,1] = paste(bp[,1], 'meta', sep='|')
bp[,2] = paste(bp[,2], 'tmp', sep='|')

g.bp = graph_from_edgelist(el = bp[,1:2], directed = FALSE)
E(g.bp)$weight = as.numeric(bp[,3])
V(g.bp)$node_type = extract_string(V(g.bp)$name, '\\|', 2)
V(g.bp)$name = extract_string(V(g.bp)$name, '\\|', 1)
V(g.bp)$native_nom = ifelse(V(g.bp)$node_type == 'meta', 'PubChemID', 'ENSPID')
V(g.bp)$PubChemID[V(g.bp)$native_nom == 'PubChemID'] = V(g.bp)$name[V(g.bp)$native_nom == 'PubChemID']
V(g.bp)$ENSPID[V(g.bp)$native_nom == 'ENSPID'] = V(g.bp)$name[V(g.bp)$native_nom == 'ENSPID']

bp.meta = vector('list', 3)
names(bp.meta) = paste('meta', names(g.nise), sep=';')
for(i in seq_along(bp.meta)){
  id = which(V(g.bp)$name %in% c(V(g.meta)$name, V(g.strings[[i]])$name, V(g.nise[[i]])$name))
  bp.meta[[i]] = induced_subgraph(graph = g.bp, vids = id)
  V(bp.meta[[i]])$node_type[V(bp.meta[[i]])$node_type != 'meta'] = names(g.nise)[i]
}

bp.ppi = vector('list', 3)
names(bp.ppi) = c('t;p', 't;ph', 'p;ph')
for(i in seq_along(bp.ppi)){
  n1 = extract_string(names(bp.ppi)[i], ';', 1)
  n2 = extract_string(names(bp.ppi)[i], ';', 2)
  v1 = unique(c(V(g.strings[[n1]])$name, V(g.nise[[n1]])$name))
  v2 = unique(c(V(g.strings[[n2]])$name, V(g.nise[[n2]])$name))
  common = intersect(v1, v2)
  
  bp.tmp = matrix(c(paste(common, n1, sep='|'), paste(common, n2, sep='|')), ncol=2)
  bp.tmp = graph_from_edgelist(el = bp.tmp, directed = FALSE)
  V(bp.tmp)$node_type = extract_string(V(bp.tmp)$name, '\\|', 2)
  V(bp.tmp)$name = extract_string(V(bp.tmp)$name, '\\|', 1)
  V(bp.tmp)$native_nom = 'ENSPID'
  V(bp.tmp)$ENSPID = V(bp.tmp)$name
  E(bp.tmp)$weight = ppi_edge_threshold
  bp.ppi[[i]] = bp.tmp
}

names(g.strings) = paste(names(g.strings), 'string', sep='_')
names(g.nise) = paste(names(g.nise), 'nise', sep='_')

input_graphs = c(list(meta = g.meta),
                 g.strings,
                 g.nise,
                 bp.meta,
                 bp.ppi)

#==================================#
# Create data list object for AMEND
#==================================#
# Names of data vectors must correspond to 'name' v.attr of graph

eci.ensp = lapply(ECI[1], function(x){
  id = match(names(x), symbol2ensp[,2])
  names(x)[!is.na(id)] = symbol2ensp[id[!is.na(id)], 1]
  x
})

eci.refseq = lapply(ECI[2:3], function(x){
  id = match(names(x), acc2ensp[,1])
  names(x)[!is.na(id)] = acc2ensp[id[!is.na(id)], 2]
  x
})

eci = c(eci.ensp, eci.refseq, ECI[4])

#======#
# AMEND ----
#======#
# AMEND functions
source(paste0(path.amend, "run_AMEND.R"))
source(paste0(path.amend, "create_integrated_graph.R"))
source(paste0(path.amend, "utils.R"))
source(paste0(path.amend, "RandomWalk.R"))
source(paste0(path.amend, "heinz.R"))

# Jump & Switch parameters
jp = c(meta = 0.5,
       t = 0.5,
       p = 0.5,
       ph = 0.5)
sp = list(t = c(string=0.5,
                nise=0.5),
          p = c(string=0.5,
                nise=0.5),
          ph = c(string=0.5,
                nise=0.5))

# Component and layer seed weights
nw = c(meta = 0.25,
       t = 0.25,
       p = 0.25,
       ph = 0.25) 
lw = list(t = c(string=0.5,
                nise=0.5),
          p = c(string=0.5,
                nise=0.5),
          ph = c(string=0.5,
                 nise=0.5))

seed.function = 'shift_scale'
fun.params = list(DOI = 1, w = 0.5)
N = 150
agg.layer = c('string', 'nise')[2]
aggregate.multiplex = list(primary = c(paste0("t_", agg.layer),
                                       paste0("p_", agg.layer),
                                       paste0("ph_", agg.layer)), agg.method = "mean")

# Create integrated and aggregated graphs
full.graph = create_integrated_graph(graph = input_graphs, data = eci, multiplex = TRUE, heterogeneous = TRUE, 
                                     FUN = seed.function, FUN.params = fun.params, lcc = TRUE)
agg.graph = create_aggregated_graph(graph = full.graph, control = aggregate.multiplex)


#=== Running AMEND ===#
db.adj.methods = list(NULL,
                      list(component = c('t_string', 'p_string', 'ph_string'), method = 'BS'),
                      list(component = c('t_string', 'p_string', 'ph_string'), method = 'IN'),
                      list(component = c('t_string', 'p_string', 'ph_string'), method = 'SDS'))

n_cores = 16

# No BRW, non-aggregated, DB adjustment
if(1){
  for(i in seq_along(db.adj.methods)){
    subnet = run_AMEND(graph = input_graphs, n = N, data = eci, brw.attr = NULL, aggregate.multiplex = NULL, degree.bias = db.adj.methods[[i]],
                       multiplex = TRUE, heterogeneous = TRUE, normalize = "degree", FUN = seed.function, FUN.params = fun.params,
                       jump.prob = jp, net.weight = nw, switch.layer.prob = sp, layer.weight = lw, verbose = TRUE, in.parallel = TRUE, n.cores = n_cores)
    
    # saveRDS(subnet, file = paste0(path.amend.res, "Oglcnac_n", N, "_db.", ifelse(is.null(db.adj.methods[[i]]), 'none', db.adj.methods[[i]]$method), ".rds"))
  }
}
# No BRW, aggregated multiplex, DB adjustment
if(1){
  for(i in seq_along(db.adj.methods)){
    subnet = run_AMEND(graph = input_graphs, n = N, data = eci, brw.attr = NULL, aggregate.multiplex = aggregate.multiplex, degree.bias = db.adj.methods[[i]],
                       multiplex = TRUE, heterogeneous = TRUE, normalize = "degree", FUN = seed.function, FUN.params = fun.params,
                       jump.prob = jp, net.weight = nw, switch.layer.prob = sp, layer.weight = lw, verbose = TRUE, in.parallel = TRUE, n.cores = n_cores)
    
    # saveRDS(subnet, file = paste0(path.amend.res, "Oglcnac_n", N, "_agg.", agg.layer, "_db.", ifelse(is.null(db.adj.methods[[i]]), 'none', db.adj.methods[[i]]$method), ".rds"))
  }
}


#===================#
# Processing Results ----
#===================#
if(0){
  db.adj.method = c('IN', 'BS', 'SDS', 'none')[1]
  agg.path = paste0('Oglcnac_n', N, '_agg.', agg.layer, '_db.', db.adj.method, '.rds')
  non.agg.path = paste0('Oglcnac_n', N, '_db.', db.adj.method, '.rds')
  
  
  subnets = list(readRDS(paste0(path.amend.res, non.agg.path)),
                 readRDS(paste0(path.amend.res, agg.path)))
  names(subnets) = c(paste0('db.', db.adj.method),
                     paste0('agg_db.', db.adj.method))
  
  # Changing NA pre.seed.values to 0
  subnets = lapply(subnets, function(x){
    V(x$module)$pre.seed.values[is.na(V(x$module)$pre.seed.values)] = 0
    x
  })
  
  # Adding Compound Names for metabolites
  meta.map = read.csv(file = paste0(path.data, 'refmet.csv'), header = TRUE)
  meta.map = rbind(meta.map, c('Hydroxide', rep(NA, 7), 961), c('Chloride', rep(NA, 7), 312), c('Phosphorylcholine', rep(NA, 7), 1014))
  meta.map$pubchem_cid = as.character(meta.map$pubchem_cid)
  rm.id = which(is.na(meta.map$pubchem_cid))
  meta.map = meta.map[-rm.id,]
  subnets = lapply(subnets, function(x) {
    id = match(V(x$module)$PubChemID, meta.map$pubchem_cid)
    V(x$module)$Compound_Name[!is.na(id)] = meta.map$refmet_name[id[!is.na(id)]]
    id = which(is.na(V(x$module)$Compound_Name) & V(x$module)$node_type == 'meta')
    V(x$module)$Compound_Name[id] = V(x$module)$PubChemID[id]
    x
  })
  id = match(V(agg.graph)$PubChemID, meta.map$pubchem_cid)
  V(agg.graph)$Compound_Name[!is.na(id)] = meta.map$refmet_name[id[!is.na(id)]]
  
  # Attempt to get better mapping between ENSPID and Symbols
  uniq.ensp = c(unlist(lapply(subnets, function(x) V(x$module)$ENSPID)), V(agg.graph)$ENSPID)
  # uniq.ensp = unlist(lapply(subnets, function(x) V(x$module)$ENSPID))
  uniq.ensp = unique(unname(uniq.ensp[!is.na(uniq.ensp)]))
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://jan2024.archive.ensembl.org")
  mart_data <- useDataset("mmusculus_gene_ensembl", mart = mart)
  ensp2symbol <- getBM(attributes = c("mgi_symbol", "ensembl_peptide_id"),
                       filters = "ensembl_peptide_id",
                       values = uniq.ensp,
                       mart = mart_data)
  ensp2symbol = ensp2symbol[apply(ensp2symbol,1,function(x) all(x != "")), ]
  
  subnets = lapply(subnets, function(x){
    id = match(V(x$module)$ENSPID, ensp2symbol[,2])
    V(x$module)$Symbol[!is.na(id)] = ensp2symbol[id[!is.na(id)],1]
    x
  })
  id = match(V(agg.graph)$ENSPID, ensp2symbol[,2])
  V(agg.graph)$Symbol[!is.na(id)] = ensp2symbol[id[!is.na(id)],1]
  
  #=== Create 'label' vertex attribute ===#
  subnets = lapply(subnets, function(x){
    V(x$module)$label = ifelse(is.na(V(x$module)$Symbol) & extract_string(V(x$module)$node_type, '_', 1) %in% c('t', 'p', 'ph'), V(x$module)$ENSPID, 
                               ifelse(is.na(V(x$module)$Symbol), V(x$module)$Compound_Name, V(x$module)$Symbol))
    V(x$module)$label.extra = paste(V(x$module)$label, V(x$module)$node_type, sep='|')
    V(x$module)$label.extra2 = ifelse(grepl('_', V(x$module)$node_type), paste(V(x$module)$label, extract_string(V(x$module)$node_type,'_',2), sep='|'), V(x$module)$label)
    x
  })
  
  #====================#
  # Choose final module
  #====================#
  subnet = subnets$agg_db.IN
  
  #==== Create node_label v.attr ====#
  V(subnet$module)$node_label = ifelse(V(subnet$module)$node_type == 'meta', V(subnet$module)$Compound_Name, V(subnet$module)$Symbol)
  V(agg.graph)$node_label = ifelse(V(agg.graph)$node_type == 'meta', V(agg.graph)$Compound_Name, V(agg.graph)$Symbol)
  
  #==================#
  # Get KEGG Pathways
  #==================#
  kegg_pathways = readRDS(file = paste0(path.data, "kegg_pathways.rds"))
  kp.comp = lapply(kegg_pathways, function(x) unname(c(x$GENE, x$COMPOUND)))
  names(kp.comp) = unlist(lapply(kegg_pathways, function(x) x$NAME))
  names(kp.comp) = gsub(' - Mus musculus \\(house mouse\\)', '', names(kp.comp))
  
  #==============#
  # ORA with KEGG ----
  #==============#
  library(fgsea)
  set.seed(839)
  uni = V(agg.graph)$node_label
  feature.names = V(subnet$module)$node_label
  ora.kegg = fora(pathways = kp.comp, genes = feature.names, universe = uni) %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::mutate(db = 'kegg',
                  overlapGenes = unlist(lapply(overlapGenes, function(x) paste(x, collapse = ", "))))
  
  #=========================#
  # ORA w/ Reactome Pathways ----
  #=========================#
  r_paths = readRDS(paste0(path.data, "m_pathways.RDS"))
  names(r_paths) = gsub('Mus musculus: ', '', names(r_paths))
  # Merging redundant pathways with "|" character
  if(1){ 
    tmp = unlist(lapply(r_paths, function(x) paste(sort(x), collapse = "|")))
    tmp2 = table(tmp)
    tmp3 = vector("list", length(tmp2))
    for(i in seq_along(tmp3)){
      tmp3[[i]] = r_paths[[match(names(tmp2)[i], tmp)]]
      nm = names(tmp)[tmp %in% names(tmp2)[i]]
      names(tmp3)[i] = nm[which.min(nchar(nm))]
    }
    r_paths = tmp3
  }
  
  uni = V(agg.graph)$node_label[V(subnet$module)$node_type != 'meta']
  gene.names = V(subnet$module)$Symbol[V(subnet$module)$node_type != 'meta']
  ora.react = fora(pathways = r_paths, genes = gene.names, universe = uni) %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::mutate(db = 'reactome',
                  overlapGenes = unlist(lapply(overlapGenes, function(x) paste(x, collapse = ", "))))
}











