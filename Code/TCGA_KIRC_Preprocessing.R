#========================#
# TCGA-KIRC Data Analysis
#========================#
#== Resources ==#
# TCGA-KIRC data download: http://firebrowse.org/?cohort=KIRC&download_dialog=true
# TCGA Barcode Description: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables
# Broad Firehose FAQ: https://broadinstitute.atlassian.net/wiki/spaces/GDAC/pages/844334036/FAQ
# Broad Firehose Documentation: https://broadinstitute.atlassian.net/wiki/spaces/GDAC/pages/844334346/Documentation
# Limma paired sample design: pg. 25 of 'design_contrast_matrices.pdf' 

#== General Info ==#
# Sample IDs:
# - Ex: TCGA-CZ-5445-01
#     TCGA = project; CZ = Tissue source site (TSS); 5445 = Patient ID; 01 = Sample type
# - Sample Type codes: 01 = Primary solid tumor; 11 = solid tissue normal

#== Questions ==#
# - For DE analysis, should I include TSS as a variable to control for potential confounding?
#   - Then again, these are paired tumor-normal samples, with both samples coming from the same site, so this may be unnecessary.
# - For mRNA/miRNA expression data, should I remove or impute NA values?
# - Snap-frozen vs. FFPE samples?
# - Mature miRNA or not? 
# - For survival analysis, should I only consider expression profiles for tumor samples? (Yes)

#== Set-up ==#
rm(list=ls())

# Path names
path.data = "/path/to/data/"

# Libraries
library(data.table); library(limma); library(survival); library(missForest); library(dplyr); library(Matrix)

# Helper functions
extract_string = function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))

# Set random seed
set.seed(351824)

#==============#
# Clinical Data ----
#==============#
clin.dat = data.table::fread(file = paste0(path.data, "All_CDEs.txt"), 
                            header = TRUE) |>
  as.data.frame()

clin.dat = data.table::transpose(clin.dat)
colnames(clin.dat) = as.character(clin.dat[1,])
clin.dat = clin.dat[-1,]
clin.dat$patient_id = toupper(clin.dat$patient_id)

# Misc 
if(0){
  View(clin.dat)
  table(table(clin.dat$patient_id))
}

#=================#
# Methylation Data ----
#=================#
if(1){
  # Mean signal values per gene
  meth.dat = data.table::fread(file = paste0(path.data, "KIRC.meth.by_mean.data.txt"),
                               header = TRUE) |>
    as.data.frame()
  
  # Misc
  if(0){
    dim(meth.dat)
    View(meth.dat)
    table(as.character(meth.dat[1,]))
    sum(is.na(meth.dat))
  }
  
  # Removing NA values
  meth.dat = na.omit(meth.dat)
  
  # Setting genes as rownames
  meth.dat = meth.dat[-1,]
  rownames(meth.dat) = meth.dat[,1]
  meth.dat = meth.dat[,-1]
  
  # Removing tissue types other than tumor or normal
  tt.tmp = substr(extract_string(colnames(meth.dat), "-", 4), 1, 2)
  meth.dat = meth.dat[,tt.tmp %in% c("01", "11")]
  
  # Getting patient IDs and tissue types
  pt.ids = extract_string(colnames(meth.dat), "-", 3)
  tt.tmp = substr(extract_string(colnames(meth.dat), "-", 4), 1, 2)
  tissue.type = ifelse(tt.tmp == "01", "tumor", ifelse(tt.tmp == "11", "normal", "other"))
  
  # Misc
  if(0){
    mean(pt.ids %in% clin.dat$patient_id)
    table(tissue.type)
    View(meth.dat)
  }
  
  # Changing meth.dat to a numeric matrix b.values
  tmp = matrix(0, nrow = nrow(meth.dat), ncol = ncol(meth.dat), dimnames = list(rownames(meth.dat), colnames(meth.dat)))
  for(i in seq_len(ncol(meth.dat))){
    tmp[,i] = as.numeric(meth.dat[,i])
  }
  b.values = tmp
  
  # Creating matrix of M-values
  m.values = log(b.values / (1 - b.values))
  
  # Comparing M and Beta values
  if(0){
    plot(density(b.values))
    plot(density(m.values))
  }
  
  # Differential methylation analysis
  # only get samples that have repeated measures to calculate corr
  multi.ids = pt.ids[pt.ids %in% names(table(pt.ids))[table(pt.ids) > 1]]
  tt.tmp = tissue.type[pt.ids %in% multi.ids]
  design <- model.matrix(~0+tt.tmp)
  dat.tmp = m.values[,pt.ids %in% multi.ids]
  # corr <- limma::duplicateCorrelation(dat.tmp, design, block = multi.ids)
  # corr = corr$consensus.correlation # 0.1103251 # Takes a long time to run
  corr = 0.1103251
  
  # Mixed effects model
  design <- model.matrix(~0+tissue.type)
  fit = lmFit(object = m.values, design = design, block = pt.ids, correlation = corr)
  
  contrast_matrix = makeContrasts("tissue.typetumor-tissue.typenormal", levels = colnames(design))
  fit = contrasts.fit(fit, contrasts = contrast_matrix)
  fit_eb = eBayes(fit, trend = TRUE)
  meth.res = topTable(fit_eb, number = 100000, adjust.method = "BH")
  
  # Misc
  if(0){
    head(meth.res)
    summary(meth.res$adj.P.Val)
    summary(meth.res$logFC)
    mean(meth.res$logFC < 0) 
  }
  
  # Visualize results
  if(0){
    genes = 1:4
    top = rownames(meth.res)[genes]
    for(i in seq_along(genes)){
      dat = b.values[rownames(b.values) == top[i],]
      g = as.factor(tissue.type)
      plot(g, dat, ylab = "Expression Values", xlab = "Groups")
    }
    meth.res[genes,]
  }
}

#=====================#
# mRNA Expression Data ----
#=====================#
if(1){
  # RSEM, quantile normalized, log2 transformed
  mrna.dat = data.table::fread(file = paste0(path.data, "KIRC.uncv2.mRNAseq_RSEM_normalized_log2.txt"),
                               header = TRUE) |>
    as.data.frame()
  
  # Misc
  if(0){
    dim(mrna.dat)
    View(mrna.dat)
    
    tmp = na.omit(mrna.dat)
    nrow(mrna.dat) - nrow(tmp)
    1-nrow(tmp)/nrow(mrna.dat)
    na.counts = apply(mrna.dat, 1, function(x) sum(is.na(x)))
    table(na.counts)
  }
  
  # Removing genes with NAs
  mrna.dat = na.omit(mrna.dat)
  
  # Remove genes without gene symbol
  mrna.dat = mrna.dat[!grepl("\\?", mrna.dat$gene),]
  
  # Gene symbol <--> Entrez gene ID mapping
  symbol_entrez.map = data.frame(Symbol = extract_string(mrna.dat$gene, "\\|", 1), Entrez = extract_string(mrna.dat$gene, "\\|", 2))
  mrna.dat$gene = extract_string(mrna.dat$gene, "\\|", 1)
  
  # Setting genes as rownames and removing gene column
  rownames(mrna.dat) = mrna.dat$gene
  mrna.dat = mrna.dat |> dplyr::select(-gene)
  
  # Removing tissue types other than tumor or normal
  mrna.dat = mrna.dat[,extract_string(colnames(mrna.dat), "-", 4) %in% c("01", "11")]
  
  # Getting patient IDs and tissue types
  pt.ids = extract_string(colnames(mrna.dat), "-", 3)
  tissue.type = ifelse(extract_string(colnames(mrna.dat), "-", 4) == "01", "tumor", ifelse(extract_string(colnames(mrna.dat), "-", 4) == "11", "normal", "other"))
  
  # Misc 
  if(0){
    mean(pt.ids %in% clin.dat$patient_id)
    table(tissue.type)
    table(table(pt.ids))
    View(mrna.dat)
  }
  
  #=== DE Analysis ===#
  # only get samples that have repeated measures to calculate corr
  multi.ids = pt.ids[pt.ids %in% names(table(pt.ids))[table(pt.ids) > 1]]
  tt.tmp = tissue.type[pt.ids %in% multi.ids]
  design <- model.matrix(~0+tt.tmp)
  dat.tmp = mrna.dat[,pt.ids %in% multi.ids]
  corr <- limma::duplicateCorrelation(dat.tmp, design, block = multi.ids)
  corr$consensus.correlation
  
  # Mixed effects model
  design <- model.matrix(~0+tissue.type)
  fit = lmFit(object = mrna.dat, design = design, block = pt.ids, correlation = corr$consensus.correlation)
  
  contrast_matrix = makeContrasts("tissue.typetumor-tissue.typenormal", levels = colnames(design))
  fit = contrasts.fit(fit, contrasts = contrast_matrix)
  fit_eb = eBayes(fit, trend = TRUE)
  mrna.res = topTable(fit_eb, number = 100000, adjust.method = "BH")
  
  # Misc
  if(0){
    head(mrna.res)
    summary(mrna.res$adj.P.Val)
    summary(mrna.res$logFC) 
  }
  
  # Visualize results
  if(0){
    genes = 1:4
    top = rownames(mrna.res)[genes]
    for(i in seq_along(genes)){
      dat = mrna.dat[rownames(mrna.dat) == top[i],]
      g = as.factor(tissue.type)
      plot(g, dat, ylab = "Expression Values", xlab = "Groups")
    }
    mrna.res[genes,]
  }
  
  # DE analysis with patient IDs as fixed effect
  if(0){
    design <- model.matrix(~0+tissue.type+pt.ids)
    fit = lmFit(object = mrna.dat, design = design)
    
    contrast_matrix = makeContrasts("tissue.typetumor-tissue.typenormal", levels = colnames(design))
    fit = contrasts.fit(fit, contrasts = contrast_matrix)
    fit_eb = eBayes(fit, trend = TRUE)
    tmp = topTable(fit_eb, number = 100000, adjust.method = "BH")
    tmp = tmp[match(rownames(mrna.res), rownames(tmp)),]
    head(tmp)
    head(mrna.res)
    summary(tmp$adj.P.Val - mrna.res$adj.P.Val)
    summary(tmp$logFC - mrna.res$logFC)
    # Not much difference. log fold changes are slightly larger and p-values are slightly lower for mixed effect model
    # So mixed effect model seems to be more powerful.
  }
}

#======================#
# miRNA Expression Data ----
#======================#
if(1){
  mature = TRUE
  if(mature){
    # RPM, log2 transformed
    mirna.dat = data.table::fread(file = paste0(path.data, "KIRC.miRseq_mature_RPM_log2.txt"),
                                  header = TRUE) |>
      as.data.frame()
  }else{
    # RPKM, log2 transformed
    mirna.dat = data.table::fread(file = paste0(path.data, "KIRC.miRseq_RPKM_log2.txt"),
                                  header = TRUE) |>
      as.data.frame()
  }
  
  # Misc
  if(0){
    # Mature: 2588 rows
    # Not Mature: 462 rows
    dim(mirna.dat)
    View(mirna.dat)
    
    tmp = na.omit(mirna.dat)
    nrow(mirna.dat) - nrow(tmp)
    1-nrow(tmp)/nrow(mirna.dat) # 90% of mature miRNAs have NAs
    na.counts = apply(mirna.dat, 1, function(x) sum(is.na(x)))
    tnac = table(na.counts)
    n.tnac = as.numeric(names(tnac))
    k = 0.75
    sum(tnac[n.tnac/ncol(mirna.dat) <= k])
  }
  
  rownames(mirna.dat) = mirna.dat[,1]
  mirna.dat = mirna.dat[,-1]
  
  # Removing tissue types other than tumor or normal
  tt.tmp = substr(extract_string(colnames(mirna.dat), "-", 4), 1, 2)
  mirna.dat = mirna.dat[,tt.tmp %in% c("01", "11")]
  
  # Getting patient IDs and tissue types
  pt.ids = extract_string(colnames(mirna.dat), "-", 3)
  tt.tmp = substr(extract_string(colnames(mirna.dat), "-", 4), 1, 2)
  tissue.type = ifelse(tt.tmp == "01", "tumor", ifelse(tt.tmp == "11", "normal", "other"))
  
  # Removing miRNAs with percent of NA samples above a certain threshold
  na.counts = apply(mirna.dat, 1, function(x) sum(is.na(x)))
  na.perc = na.counts / ncol(mirna.dat)
  k = 0.25
  # sum(na.perc <= k)
  mirna.dat = mirna.dat[na.perc <= k,]
  
  # Random forest imputation for remaining missing values
  impute = TRUE
  if(impute){
    mirna.dat.t = t(mirna.dat)
    mirna.dat.t.imp = missForest(mirna.dat.t)
    mirna.dat.imp = t(mirna.dat.t.imp$ximp)
    mirna.dat0 = mirna.dat.imp
  }else mirna.dat0 = mirna.dat
  
  #=== DE Analysis ===#
  # only get samples that have repeated measures to calculate corr
  multi.ids = pt.ids[pt.ids %in% names(table(pt.ids))[table(pt.ids) > 1]]
  tt.tmp = tissue.type[pt.ids %in% multi.ids]
  design <- model.matrix(~0+tt.tmp)
  dat.tmp = mirna.dat0[,pt.ids %in% multi.ids]
  corr <- limma::duplicateCorrelation(dat.tmp, design, block = multi.ids)
  corr$consensus.correlation
  
  # Mixed effects model
  design <- model.matrix(~0+tissue.type)
  fit = lmFit(object = mirna.dat0, design = design, block = pt.ids, correlation = corr$consensus.correlation)
  
  contrast_matrix = makeContrasts("tissue.typetumor-tissue.typenormal", levels = colnames(design))
  fit = contrasts.fit(fit, contrasts = contrast_matrix)
  fit_eb = eBayes(fit, trend = TRUE) # Assuming library sizes of samples are fairly similar
  mirna.res = topTable(fit_eb, number = 100000, adjust.method = "BH")
  
  # Misc
  if(0){
    head(mirna.res)
    summary(mirna.res$adj.P.Val)
    summary(mirna.res$logFC)
  }
  
  # Visualize results
  if(0){
    genes = 1:4
    top = rownames(mirna.res)[genes]
    for(i in seq_along(genes)){
      dat = mirna.dat[rownames(mirna.dat) == top[i],]
      g = as.factor(tissue.type)
      plot(g, dat, ylab = "Expression Values", xlab = "Groups")
    }
    mirna.res[genes,]
  }
  
  # DE analysis with non-imputed data
  if(0){
    mirna.dat.narm = na.omit(mirna.dat)
    
    # only get samples that have repeated measures to calculate corr
    multi.ids = pt.ids[pt.ids %in% names(table(pt.ids))[table(pt.ids) > 1]]
    tt.tmp = tissue.type[pt.ids %in% multi.ids]
    design <- model.matrix(~0+tt.tmp)
    dat.tmp = mirna.dat.narm[,pt.ids %in% multi.ids]
    corr <- limma::duplicateCorrelation(dat.tmp, design, block = multi.ids)
    corr$consensus.correlation
    
    # Mixed effects model
    design <- model.matrix(~0+tissue.type)
    fit = lmFit(object = mirna.dat.narm, design = design, block = pt.ids, correlation = corr$consensus.correlation)
    
    contrast_matrix = makeContrasts("tissue.typetumor-tissue.typenormal", levels = colnames(design))
    fit = contrasts.fit(fit, contrasts = contrast_matrix)
    fit_eb = eBayes(fit, trend = TRUE) # Assuming library sizes of samples are fairly similar
    mirna.res.narm = topTable(fit_eb, number = 100000, adjust.method = "BH")
    ids = match(rownames(mirna.res), rownames(mirna.res.narm))
    ids = ids[!is.na(ids)]
    mirna.res.narm = mirna.res.narm[ids,]
    tmp = mirna.res[rownames(mirna.res) %in% rownames(mirna.res.narm),]
    head(mirna.res.narm)
    head(tmp)
    summary(mirna.res.narm$adj.P.Val - tmp$adj.P.Val)
    summary(mirna.res.narm$logFC - tmp$logFC)
    # Very little difference 
    
    tmp = mirna.res[!rownames(mirna.res) %in% rownames(mirna.res.narm),]
    summary(tmp$logFC)
    summary(tmp$adj.P.Val)
    mean(tmp$adj.P.Val <= 0.05)
  }
  
  # DE analysis with patient IDs as fixed effect
  if(0){
    design <- model.matrix(~0+tissue.type+pt.ids)
    fit = lmFit(object = mirna.dat0, design = design)
    
    contrast_matrix = makeContrasts("tissue.typetumor-tissue.typenormal", levels = colnames(design))
    fit = contrasts.fit(fit, contrasts = contrast_matrix)
    fit_eb = eBayes(fit, trend = TRUE)
    tmp = topTable(fit_eb, number = 100000, adjust.method = "BH")
    tmp = tmp[match(rownames(mirna.res), rownames(tmp)),]
    head(tmp)
    head(mirna.res)
    summary(tmp$adj.P.Val - mirna.res$adj.P.Val)
    summary(tmp$logFC - mirna.res$logFC)
    # Not much difference. log fold changes are slightly larger and p-values are slightly lower for fixed effect model
    # So fixed effect model seems to be more powerful...
  }
}

#=======#
# Saving
#=======#
if(0){
  # Putting Omics data in one object
  omics.dat = list(mrna = mrna.dat, mirna = mirna.dat0, methyl = m.values)
  saveRDS(omics.dat, file = "/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Results/multi_omics.rds")
  
  # Putting DE results in one object
  de.res = list(mrna = mrna.res, mirna = mirna.res, methyl = meth.res)
  saveRDS(de.res, file = "/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Results/DE_results.rds")
}

# omics.dat = readRDS("/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Results/multi_omics.rds")
# de.res = readRDS("/Users/samboyd/Documents/GRA/Network Analysis/Paper 3/Results/DE_results.rds")

if(0){
  lapply(omics.dat, function(x){
    c.names = colnames(x)
    tmp = table(unlist(lapply(strsplit(c.names, '-'), function(y) y[3])))
    table(tmp)
  })
}

#==================#
# Survival Analysis ----
#==================#
surv.dat = clin.dat[,c("patient_id", "vital_status", "days_to_death")]
# Change vital_status to 0=alive, 1=dead
surv.dat$vital_status = ifelse(surv.dat$vital_status == "alive", 0, 1)
surv.dat$days_to_death = as.numeric(surv.dat$days_to_death)
# Getting only tumor samples
omics.dat = lapply(omics.dat, function(x){
  tt.tmp = substr(extract_string(colnames(x), "-", 4), 1, 2)
  x[,tt.tmp == "01"]
})
# Extracting omics data and appending to surv.dat
surv.dat.list = vector("list", length(omics.dat)); names(surv.dat.list) = names(omics.dat)
r.ids = surv.dat$patient_id
for(i in seq_along(omics.dat)){
  message(i)
  omics.res = matrix(0, nrow = length(r.ids), ncol = nrow(omics.dat[[i]]), dimnames = list(r.ids, rownames(omics.dat[[i]])))
  pt.ids = extract_string(colnames(omics.dat[[i]]), "-", 3)
  for(j in seq_along(r.ids)){
    if(r.ids[j] %in% pt.ids){
      omics.res[j,] = omics.dat[[i]][,pt.ids == r.ids[j]]
    }else{
      omics.res[j,] = rep(NA, ncol(omics.res))
    }
  }
  # Get rows in same order as surv.dat, then cbind with surv.dat
  omics.res = omics.res[match(r.ids, rownames(omics.res)),]
  surv.dat.list[[i]] = cbind(surv.dat, omics.res)
  surv.dat.list[[i]] = surv.dat.list[[i]][!is.na(omics.res[,1]),]
  # Scaling by standard deviation feature-wise before cox PH modeling (no need to center since survival::coxph does this for you)
  omics.res = as.matrix(surv.dat.list[[i]][,-c(1:ncol(surv.dat))])
  surv.dat.list[[i]][,-c(1:ncol(surv.dat))] = omics.res %*% Matrix::Diagonal(x = apply(omics.res, 2, sd)^(-1))
}

# Misc
if(0){
  lapply(surv.dat.list, function(x) x[1:10,1:10])
  lapply(surv.dat.list, dim)
  lapply(omics.dat, dim)
  
  # Checking for special characters that may interfere with formula call
  lapply(surv.dat.list, function(x){
    v = colnames(x)[4:ncol(x)]
    id = grep(pattern = "-|\\+", x = v)
    v[id]
  })
}

# Modify feature names to not conflict with formula in coxph
feature.names = vector("list", length(surv.dat.list))
for(i in seq_along(feature.names)){
  res = matrix("", nrow = ncol(surv.dat.list[[i]])-3, ncol = 2, dimnames = list(1:(ncol(surv.dat.list[[i]])-3), c("Original", "Modified")))
  v = colnames(surv.dat.list[[i]])[4:ncol(surv.dat.list[[i]])]
  res[,1] = v
  res[,2] = make.names(v, unique = TRUE)
  feature.names[[i]] = res
  colnames(surv.dat.list[[i]])[4:ncol(surv.dat.list[[i]])] = res[,2]
}

# For each omic type, for each feature, fit a separate Cox PH model. vital_status=1 represents a death.
# These Cox PH models will give larger coefficient estimates to features whose increase in expression is associated with an increase hazard of death.
cph.models = vector("list", length(surv.dat.list))
names(cph.models) = names(surv.dat.list)
for(i in seq_along(surv.dat.list)){
  cph.res = vector("list", ncol(surv.dat.list[[i]])-3)
  names(cph.res) = colnames(surv.dat.list[[i]])[4:ncol(surv.dat.list[[i]])]
  message(length(cph.res))
  for(j in seq_along(cph.res)){
    if(j %% 500 == 0) message(paste0("i",i,".j",j))
    v = colnames(surv.dat.list[[i]])[j+3]
    form = as.formula(paste("Surv(days_to_death, vital_status)~", v))
    cph.res[[j]] = coxph(form, data = surv.dat.list[[i]])
  }
  cph.models[[i]] = cph.res
}

# Replace correct feature names
for(i in seq_along(surv.dat.list)){
  colnames(surv.dat.list[[i]])[4:ncol(surv.dat.list[[i]])] = feature.names[[i]][, "Original"]
  names(cph.models[[i]]) = feature.names[[i]][, "Original"]
}
# saveRDS(cph.models, file = "cph_models.rds")
# saveRDS(surv.dat.list, file = 'surv_expr_data.rds')

# Misc
if(0){
  cph.models[[1]][[1]]
  summary(cph.models[[1]][[1]])$coefficients
}

# Extract hazard ratios and p-values from coxph models
hr = lapply(cph.models, function(x) lapply(x, function(y) {
  tmp = summary(y)$coefficients[,c(2,5)]
  names(tmp) = c("hr", "pval")
  tmp
  }))

# saveRDS(hr, file = "hazard_ratios.rds")
