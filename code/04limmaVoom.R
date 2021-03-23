###########################################################################
### Study gene expression upon PTEN status in TCGA
### Eddie Imada & Luigi Marchionni
### Collaboration with Tamara Lotan
### Linear model analysis
###########################################################################
### Set wd
setwd("~/Dropbox (MechPred)/Projects/PtenERG/Manuscript/")

### Clean
rm(list=ls())

### Load libraries
library(Biobase)
library(limma)
require(edgeR)
require(sva)
library(dplyr)
library(RColorBrewer)
library(gplots)
###########################################################################
### load
load("objs/pradHighGistic.rda")
######################################################################
### Create groups
pten.stat <- factor(prad$PTENstatusCNV)
erg.stat <- factor(prad$ERGstatus)
groups <- paste(pten.stat, erg.stat, sep="_")
groups <- gsub(" ", "_", groups)
groups <- factor(groups)
### Model matrix without intectept
dMat <- model.matrix( ~0+groups)
colnames(dMat) <- gsub("groups", "", colnames(dMat))
colnames(dMat) <- gsub(" ", "_", colnames(dMat))

### Make contrast matrix
cMat <- makeContrasts(
    levels=colnames(dMat),
    ## PTEN
    PTEN_NEGvsPTEN_POS= ((PTEN_NEGATIVE_ERG_NEGATIVE - PTEN_POSITIVE_ERG_NEGATIVE) + (PTEN_NEGATIVE_ERG_POSITIVE - PTEN_POSITIVE_ERG_POSITIVE))/2,
    PTEN_NEGvsPTEN_POSinERGpos = PTEN_NEGATIVE_ERG_POSITIVE - PTEN_POSITIVE_ERG_POSITIVE,
    PTEN_NEGvsPTEN_POSinERGneg = PTEN_NEGATIVE_ERG_NEGATIVE - PTEN_POSITIVE_ERG_NEGATIVE,
    ERG_POSvsERG_NEG = ((PTEN_NEGATIVE_ERG_POSITIVE + PTEN_POSITIVE_ERG_POSITIVE)/2) - ((PTEN_NEGATIVE_ERG_NEGATIVE + PTEN_POSITIVE_ERG_NEGATIVE)/2)
    )

###########################################################################
###########################################################################
### GLM using voom() and limma rather than the GLM approach in edgeR

### Prepare expression
exp <- assays(prad)[[1]]
dge <- as.matrix(round(exp))

### Create DGEList
dge <- DGEList(counts=dge, group=groups, genes=rowData(prad))
dge <- calcNormFactors(dge, method="TMMwsp")
### Filter low counts by expression min CPM 0.1 across all samples of lowest level & min total count = 15
keepGns <- filterByExpr(dge, group = pten.stat, min.count = 5)
table(keepGns)
dge <- dge[keepGns, ]
de <- voom(dge, design = dMat, normalize.method = "none",
           plot = TRUE)

######################################################################
### Without sva, using voomed data and limma: Fit glm
fit <- lmFit(de, dMat)
fit <- contrasts.fit(fit, cMat)
fit <- eBayes(fit)

### Check number of DE genes
summary(decideTests(fit, p.value = 0.01))

###########################################################################
### Extract results
tG <- lapply(colnames(fit$coefficients), function(x, y) {
    topTable(y, coef=x, n=Inf, genelist = fit$genes, resort.by="t")
}, y=fit)
names(tG) <- colnames(fit$coefficients)

#### Since genes are not redundant
tGnr <- tG
###########################################################################
### Save
cutoff <- 0.01
### Select columns
keep <- c("geneID", "geneName", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "geneType", "CAT_geneClass", "CAT_DHS_type", "HGNC_symbol", "HGNC_name", "entrez_ID", "lncrnadb_ID", "GENCODE_ID", "GENCODE_SYMBOL")
### Write DGE tables to csv files
lapply(1:length(tGnr), function(i) {
    tGnr[[i]] <- filter(tGnr[[i]], adj.P.Val < cutoff)
    tGnr[[i]] <- arrange(tGnr[[i]], desc(abs(logFC)))
    write.csv(tGnr[[i]][,keep],
              file = paste0("text/", "DGE_final_", names(tGnr[i]), ".csv"),
              row.names = FALSE)
})

### Save objs
save(tGnr, file="objs/DGE_prad.rda")

###########################################################################
### Session information and clean quit

### Session
sessionInfo()

### Clean
rm(list=ls())

