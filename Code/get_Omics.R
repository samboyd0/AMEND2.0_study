#============================================#
# Omics Data for Degree Bias and EHR Analysis
#============================================#
# 5 gene expression datasets 
# GSE112680: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112680
# GSE30219: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30219
# GSE75214 (2): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75214
# GSE3790: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3790

# limma: https://www.bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
# page 70 for eBayes function
# robust=TRUE: https://support.bioconductor.org/p/56560/#:~:text=The%20major%20purpose%20of%20arrayWeights,to%20remove%20from%20the%20analysis.

# Path names
path.data = "/path/to/data/"
path.results = "/path/to/results/"

# DE analysis
library(GEOquery)
library(biomaRt)
library(limma)
library(stringr)

#==========#
# GSE112680 ----
#==========#
# Platform: GPL10558, Illumina HumanHT-12 V4.0 expression beadchip
# Expression values have been log2 transformed and quantile normalized
als <- getGEO(GEO = "GSE112680", GSEMatrix = TRUE, getGPL = FALSE)
als = als[[1]]

pd = pData(als); dim(pd)
ex = exprs(als); dim(ex)
# sum(str_detect(rownames(ex), "NA"))

fd = data.table::fread(paste0(path.data, "GPL10558-50081.txt")); dim(fd)
mean(fd$ID %in% rownames(ex))
fd = fd[fd$ID %in% rownames(ex),]
mean(fd$Entrez_Gene_ID %in% c("", NA))
length(unique(fd$Entrez_Gene_ID[!fd$Entrez_Gene_ID %in% c("", NA)]))

# Removing Mimick samples
keep.acc = pd$geo_accession[pd$source_name_ch1 != "whole blood_MIM"]
keep.pd.idx = which(pd$geo_accession %in% keep.acc)
keep.ex.idx = which(colnames(ex) %in% keep.acc)

pd = pd[keep.pd.idx,]
ex = ex[,keep.ex.idx]

pd$group = ifelse(pd$source_name_ch1 == "whole blood_CON", "CTRL", "ALS")

# Create design matrix
X = model.matrix(~0+pd$group)
colnames(X) = c("ALS", "CTRL")

# Create contrast matrix
cnt_mat = makeContrasts(contrasts = "ALS-CTRL", levels = colnames(X))

fit = lmFit(ex, design = X)
fit = contrasts.fit(fit, contrasts = cnt_mat)
fit = eBayes(fit, trend = TRUE, robust = TRUE)
plotSA(fit) # This justifies use of trend=TRUE
res.als = topTable(fit, number = Inf, adjust.method = "BH")

# Mapping probes to Entrez IDs
logfc = res.als$logFC
names(logfc) = fd$Entrez_Gene_ID[match(rownames(res.als), fd$ID)]
logfc = logfc[!names(logfc) %in% c("", NA)]
# sum(do.call("c", lapply(strsplit(names(logfc), "///"), function(x) length(x) > 1))) # All probes map to a single Entrez ID
# length(unique(names(logfc))) != length(logfc) # Multiple probes map to the same Entrez ID
logfc = aggregate(x = logfc, by = list(names(logfc)), FUN = median)

pval = res.als$P.Value
names(pval) = fd$Entrez_Gene_ID[match(rownames(res.als), fd$ID)]
pval = pval[!names(pval) %in% c("", NA)]
pval = aggregate(x = pval, by = list(names(pval)), FUN = max)

dat.als = data.frame(logFC = logfc$x, P.Value = pval$x, row.names = logfc$Group.1)

saveRDS(dat.als, file = paste0(path.results, "ALS_DE_results.rds"))

#==========#
# GSE30219 ----
#==========#
# Platform: GPL570, [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# Only consider non-small cell lung cancer
# Data was normalized by Robust Multi-Array average (background-corrected, normalized, and log2 transformed)
lc <- getGEO("GSE30219", GSEMatrix = TRUE, getGPL = TRUE) # May not work first time with getGPL=TRUE
if (length(lc) > 1) idx <- grep("GPL570", attr(lc, "names")) else idx <- 1
lc <- lc[[idx]]

pd = pData(lc)
ex = exprs(lc)
fd = fData(lc)
length(unique(fd$ENTREZ_GENE_ID[!fd$ENTREZ_GENE_ID %in% c("", NA)]))

# Removing Small-cell samples
keep.acc = pd$geo_accession[!pd$characteristics_ch1.3 %in% c("histology: Other", "histology: SCC")]
keep.pd.idx = which(pd$geo_accession %in% keep.acc)
keep.ex.idx = which(colnames(ex) %in% keep.acc)

pd = pd[keep.pd.idx,]
ex = ex[,keep.ex.idx]

pd$group = ifelse(pd$characteristics_ch1.3 == "histology: NTL", "CTRL", "LC")

# Create design matrix
X = model.matrix(~0+pd$group)
colnames(X) = c("CTRL", "LC")

# Create contrast matrix
cnt_mat = makeContrasts(contrasts = "LC-CTRL", levels = colnames(X))

fit = lmFit(ex, design = X)
fit = contrasts.fit(fit, contrasts = cnt_mat)
fit = eBayes(fit, trend = TRUE, robust = TRUE)
# plotSA(fit) # This justifies use of trend=TRUE
res.lc = topTable(fit, number = Inf, adjust.method = "BH")

# Mapping probes to Entrez IDs
logfc = res.lc$logFC
names(logfc) = fd$ENTREZ_GENE_ID[match(rownames(res.lc), fd$ID)]
logfc = logfc[!names(logfc) %in% c("", NA)]
sum(do.call("c", lapply(strsplit(names(logfc), "///"), function(x) length(x) > 1)))
logfc = logfc[do.call("c", lapply(strsplit(names(logfc), "///"), function(x) length(x) == 1))]
length(unique(names(logfc))) != length(logfc) # Multiple probes map to the same Entrez ID
logfc = aggregate(x = logfc, by = list(names(logfc)), FUN = median)

pval = res.lc$P.Value
names(pval) = fd$ENTREZ_GENE_ID[match(rownames(res.lc), fd$ID)]
pval = pval[!names(pval) %in% c("", NA)]
pval = pval[do.call("c", lapply(strsplit(names(pval), "///"), function(x) length(x) == 1))]
pval = aggregate(x = pval, by = list(names(pval)), FUN = max)

dat.lc = data.frame(logFC = logfc$x, P.Value = pval$x, row.names = logfc$Group.1)

saveRDS(dat.lc, file = paste0(path.results, "LC_DE_results.rds"))

#==========#
# GSE75214 (2) ----
#==========#
# Platform: GPL6244	[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]
# Data was normalized by Robust Multi-Array average (background-corrected, normalized, and log2 transformed)
ibd <- getGEO("GSE75214", GSEMatrix = TRUE, getGPL = TRUE, AnnotGPL = TRUE)
if (length(ibd) > 1) idx <- grep("GPL6244", attr(ibd, "names")) else idx <- 1
ibd <- ibd[[idx]]

pd = pData(ibd)
ex = exprs(ibd)
fd = fData(ibd)
length(unique(do.call("c", strsplit(fd$`Gene ID`[fd$`Gene ID` != ""], "///"))))

# Separating into Ulcerative Colitis (UC)
uc.acc = pd$geo_accession[!pd$characteristics_ch1.1 %in% c("disease: CD", "disease: Crohn's disease")]
uc.pd.idx = which(pd$geo_accession %in% uc.acc)
uc.ex.idx = which(colnames(ex) %in% uc.acc)
uc.pd = pd[uc.pd.idx,]
uc.ex = ex[,uc.ex.idx]
uc.pd$group = ifelse(uc.pd$characteristics_ch1.1 == "disease: control", "CTRL", "UC")

# Create design matrix
X = model.matrix(~0+uc.pd$group)
colnames(X) = c("CTRL", "UC")

# Create contrast matrix
cnt_mat = makeContrasts(contrasts = "UC-CTRL", levels = colnames(X))

fit = lmFit(uc.ex, design = X)
fit = contrasts.fit(fit, contrasts = cnt_mat)
fit = eBayes(fit, trend = TRUE, robust = FALSE)
# plotSA(fit)
res.uc = topTable(fit, number = Inf, adjust.method = "BH")

# Mapping probes to Entrez IDs
logfc = res.uc$logFC
names(logfc) = fd$`Gene ID`[match(rownames(res.uc), fd$ID)]
logfc = logfc[!names(logfc) %in% c("", NA)]
sum(do.call("c", lapply(strsplit(names(logfc), "///"), function(x) length(x) > 1))) # All probes map to a single Entrez ID
logfc = logfc[do.call("c", lapply(strsplit(names(logfc), "///"), function(x) length(x) == 1))]
# length(unique(names(logfc))) != length(logfc)) # Multiple probes map to the same Entrez ID
logfc = aggregate(x = logfc, by = list(names(logfc)), FUN = median)

pval = res.uc$P.Value
names(pval) = fd$`Gene ID`[match(rownames(res.uc), fd$ID)]
pval = pval[!names(pval) %in% c("", NA)]
pval = pval[do.call("c", lapply(strsplit(names(pval), "///"), function(x) length(x) == 1))]
pval = aggregate(x = pval, by = list(names(pval)), FUN = max)

dat.uc = data.frame(logFC = logfc$x, P.Value = pval$x, row.names = logfc$Group.1)

saveRDS(dat.uc, file = paste0(path.results, "UC_DE_results.rds"))

# Separating into Chron's Disease (CD)
cd.acc = pd$geo_accession[!pd$characteristics_ch1.1 %in% c("disease: ulcerative colitis")]
cd.pd.idx = which(pd$geo_accession %in% cd.acc)
cd.ex.idx = which(colnames(ex) %in% cd.acc)
cd.pd = pd[cd.pd.idx,]
cd.ex = ex[,cd.ex.idx]
cd.pd$group = ifelse(cd.pd$characteristics_ch1.1 == "disease: control", "CTRL", "CD")

# Create design matrix
X = model.matrix(~0+cd.pd$group)
colnames(X) = c("CD", "CTRL")

# Create contrast matrix
cnt_mat = makeContrasts(contrasts = "CD-CTRL", levels = colnames(X))

fit = lmFit(cd.ex, design = X)
fit = contrasts.fit(fit, contrasts = cnt_mat)
fit = eBayes(fit, trend = TRUE, robust = FALSE)
# plotSA(fit)
res.cd = topTable(fit, number = Inf, adjust.method = "BH")

# Mapping probes to Entrez IDs
logfc = res.cd$logFC
names(logfc) = fd$`Gene ID`[match(rownames(res.cd), fd$ID)]
logfc = logfc[!names(logfc) %in% c("", NA)]
sum(do.call("c", lapply(strsplit(names(logfc), "///"), function(x) length(x) > 1))) # All probes map to a single Entrez ID
logfc = logfc[do.call("c", lapply(strsplit(names(logfc), "///"), function(x) length(x) == 1))]
# length(unique(names(logfc))) != length(logfc)) # Multiple probes map to the same Entrez ID
logfc = aggregate(x = logfc, by = list(names(logfc)), FUN = median)

pval = res.cd$P.Value
names(pval) = fd$`Gene ID`[match(rownames(res.cd), fd$ID)]
pval = pval[!names(pval) %in% c("", NA)]
pval = pval[do.call("c", lapply(strsplit(names(pval), "///"), function(x) length(x) == 1))]
pval = aggregate(x = pval, by = list(names(pval)), FUN = max)

dat.cd = data.frame(logFC = logfc$x, P.Value = pval$x, row.names = logfc$Group.1)

saveRDS(dat.cd, file = paste0(path.results, "CD_DE_results.rds"))

#==========#
# GSE3790 ----
#==========#
# Platform: GPL96 & GPL97 [HG-U133A] Affymetrix Human Genome U133 A & B Arrays
# Only using array A
# MAS5-calculated signal intensities
# Values will be log-transformed in this code
# RMA vs. MAS5: https://www.biostars.org/p/7687/#7689
hd <- getGEO("GSE3790", GSEMatrix = TRUE, AnnotGPL = TRUE)
hd <- hd[[1]]

pd = pData(hd)
ex = exprs(hd)
fd = fData(hd)
length(unique(fd$`Gene ID`[!fd$`Gene ID` %in% c("", NA)]))
sum(do.call("c", lapply(strsplit(fd$`Gene ID`, "///"), function(x) length(x) > 1)))

# code from GEO2R for log-transforming data
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if(LogC){
  ex[which(ex <= 0)] <- NaN
  exprs(hd) <- log2(ex)
  ex = exprs(hd)
}

# Box plots... samples look normalized
# par(mar=c(7,4,2,1))
# boxplot(ex, boxwex = 0.6, notch = T, main = paste ("GSE3790", "/", annotation(hd), sep = ""), outline = FALSE, las = 2)

# Removing non-caudate nucleus and HD grade 0-1 samples
case.idx = which(str_detect(pd$characteristics_ch1, "human caudate nucleus HD grade 2") | str_detect(pd$characteristics_ch1, "human caudate nucleus HD grade 3") | str_detect(pd$characteristics_ch1, "human caudate nucleus HD grade 4"))
ctrl.idx = which(str_detect(pd$characteristics_ch1, "human caudate nucleus control"))

keep.acc = pd$geo_accession[c(case.idx, ctrl.idx)]
keep.pd.idx = which(pd$geo_accession %in% keep.acc)
keep.ex.idx = which(colnames(ex) %in% keep.acc)

pd = pd[keep.pd.idx,]
ex = ex[,keep.ex.idx]

pd$group = ifelse(1:nrow(pd) %in% ctrl.idx, "CTRL", "HD")

# Create design matrix
X = model.matrix(~0+pd$group)
colnames(X) = c("CTRL", "HD")

# Create contrast matrix
cnt_mat = makeContrasts(contrasts = "HD-CTRL", levels = colnames(X))

fit = lmFit(ex, design = X)
fit = contrasts.fit(fit, contrasts = cnt_mat)
fit = eBayes(fit, trend = TRUE, robust = TRUE)
# plotSA(fit)
res.hd = topTable(fit, number = Inf, adjust.method = "BH")

# Mapping probes to Entrez IDs
logfc = res.hd$logFC
names(logfc) = fd$`Gene ID`[match(rownames(res.hd), fd$ID)]
logfc = logfc[!names(logfc) %in% c("", NA)]
sum(do.call("c", lapply(strsplit(names(logfc), "///"), function(x) length(x) > 1))) # All probes map to a single Entrez ID
logfc = logfc[do.call("c", lapply(strsplit(names(logfc), "///"), function(x) length(x) == 1))]
# length(unique(names(logfc))) != length(logfc)) # Multiple probes map to the same Entrez ID
logfc = aggregate(x = logfc, by = list(names(logfc)), FUN = median)

pval = res.hd$P.Value
names(pval) = fd$`Gene ID`[match(rownames(res.hd), fd$ID)]
pval = pval[!names(pval) %in% c("", NA)]
pval = pval[do.call("c", lapply(strsplit(names(pval), "///"), function(x) length(x) == 1))]
pval = aggregate(x = pval, by = list(names(pval)), FUN = max)

dat.hd = data.frame(logFC = logfc$x, P.Value = pval$x, row.names = logfc$Group.1)

saveRDS(dat.hd, file = paste0(path.results, "HD_DE_results.rds"))







