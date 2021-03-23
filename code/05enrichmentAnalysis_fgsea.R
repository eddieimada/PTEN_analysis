###########################################################################
### Study gene expression upon PTEN status in TCGA
### Eddie Imada & Luigi Marchionni
### Collaboration with Tamara Lotan
### Linear model analysis
###########################################################################
### Load libraries
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
register(SerialParam())

# Load DGE results
load("/Users/elimada/Dropbox (MechPred)/Projects/PtenERG/Manuscript/objs/DGE_prad_all_new.rda")

# Create named vectors with rank
XDE_overall <- read.csv("/Users/elimada/Dropbox (MechPred)/Projects/PtenERG/XDE/text/PTENposVSneg.csv", stringsAsFactors = F)
XDE_overall <- setNames(XDE_overall$BES_SUM, XDE_overall$X)
XDE_overall <- XDE_overall * -1
XDE_ERGpos <- read.csv("/Users/elimada/Dropbox (MechPred)/Projects/PtenERG/XDE/text/PTENinERGPosxde.csv", stringsAsFactors = F)
XDE_ERGpos <- setNames(XDE_ERGpos$BES_SUM, XDE_ERGpos$X)
XDE_ERGpos <- XDE_ERGpos * -1
XDE_ERGneg <- read.csv("/Users/elimada/Dropbox (MechPred)/Projects/PtenERG/XDE/text/PTENinERGnegxde.csv", stringsAsFactors = F)
XDE_ERGneg <- setNames(XDE_ERGneg$BES_SUM, XDE_ERGneg$X)
XDE_ERGneg <- XDE_ERGneg * -1

TCGA_overall <- setNames(tGnr$PTEN_NEGvsPTEN_POS$t, tGnr$PTEN_NEGvsPTEN_POS$HGNC_symbol)
TCGA_overall <- TCGA_overall[!is.na(names(TCGA_overall))]

TCGA_ERGpos <- setNames(tGnr$PTEN_NEGvsPTEN_POSinERGpos$t, tGnr$PTEN_NEGvsPTEN_POSinERGpos$HGNC_symbol)
TCGA_ERGpos <- TCGA_ERGpos[!is.na(names(TCGA_ERGpos))]

TCGA_ERGneg <- setNames(tGnr$PTEN_NEGvsPTEN_POSinERGneg$t, tGnr$PTEN_NEGvsPTEN_POSinERGneg$HGNC_symbol)
TCGA_ERGneg <- TCGA_ERGneg[!is.na(names(TCGA_ERGneg))]


var <- c("XDE_overall", "XDE_ERGpos", "XDE_ERGneg", "TCGA_overall", "TCGA_ERGpos", "TCGA_ERGneg")

set.seed(42)

# Perform Enrichment Analysis
AR <- gmtPathways("/Users/elimada/Dropbox (MechPred)/Projects/PtenERG/TCGA/misc/FGS_MOUSE_ANDROGEN_RESPONSE.gmt")
sapply(var, function(x){
    res <- fgseaMultilevel(AR, get(x), nproc = 8, minSize = 10, maxSize = 1500, sampleSize = 1001,)
    res <- res[order(padj),]
    res <- res[res$padj <= 0.1]
    fwrite(res, file = paste0("text/fgsea_AR_", x, ".csv"), sep = ",")
})


biocarta <- gmtPathways("~/Dropbox (MechPred)/Databases/MSigDB/c2.cp.biocarta.v7.0.symbols.gmt")
biocarta <- biocarta[order(unlist(lapply(biocarta, length)), decreasing = T)]
l <- unlist(lapply(biocarta, length))
biocarta <- biocarta[l >= 10 & l <= 1500]
sapply(var, function(x){
    res <- fgseaMultilevel(biocarta, get(x), nproc = 8, minSize = 10, maxSize = 1500, sampleSize = 1001,)
    res <- res[order(padj),]
    print(paste0("Processing ", x, " for Biocarta"))
    keep <- collapsePathways(res, 
                             biocarta, get(x),
                             nperm = 1000,
                             pval.threshold = 0.2)
    res <- res[pathway %in% keep$mainPathways,]
    res$padj <- p.adjust(res$pval, "BH")
    res <- res[order(padj),]
    res <- res[res$padj <= 0.1]
    fwrite(res, file = paste0("text/fgsea_biocarta_", x, ".csv"), sep = ",")
})

reactome <- gmtPathways("~/Dropbox (MechPred)/Databases/MSigDB/c2.cp.reactome.v7.0.symbols.gmt")
reactome <- reactome[order(unlist(lapply(reactome, length)), decreasing = T)]
l <- unlist(lapply(reactome, length))
reactome <- reactome[l >= 10 & l <= 1500]
sapply(var, function(x){
    res <- fgseaMultilevel(reactome, get(x), nproc = 8, minSize = 10, maxSize = 1500, sampleSize = 201)
    res <- res[order(padj),]
    print(paste0("Processing ", x, " for Reactome"))
    keep <- collapsePathways(res[padj <= 0.1],
                             reactome, get(x),
                             nperm = 500,
                             pval.threshold = 0.05)
    res <- res[pathway %in% keep$mainPathways,]
    fwrite(res, file = paste0("text/fgsea_reactome_", x, ".csv"), sep = ",")
})

#############################################################################

hallmarks <- gmtPathways("~/Dropbox (MechPred)/Databases/MSigDB/h.all.v7.0.symbols.gmt")
hallmarks <- hallmarks[order(unlist(lapply(hallmarks, length)), decreasing = T)]
l <- unlist(lapply(hallmarks, length))
hallmarks <- hallmarks[l >= 15 | l <= 1500]
sapply(var, function(x){
    print(paste0("Processing ", x, " for Hallmarks"))
    res <- fgseaMultilevel(hallmarks, get(x), nproc = 8, minSize = 15, maxSize = 1500, sampleSize = 1001)
    res <- res[order(padj),]
    res <- res[res$padj <= 0.1]
    fwrite(res, file = paste0("text/fgsea_hallmarks_", x, ".csv"), sep = ",")
})

#############################################################################

BP <- gmtPathways("~/Dropbox (MechPred)/Databases/MSigDB/c5.bp.v7.0.symbols.gmt")
BP <- BP[order(unlist(lapply(BP, length)), decreasing = F)]
l <- unlist(lapply(BP, length))
BP <- BP[l >= 15 | l <= 300]
sapply(var, function(x){
    print(paste0("Processing ", x, " for BP"))
    res <- fgseaMultilevel(BP, get(x), nproc = 8, minSize = 15, maxSize = 300, sampleSize = 201)
    res <- res[order(padj),]
    keep <- collapsePathways(res[padj <= 0.1],
                             BP, get(x),
                             nperm = 200,
                             pval.threshold = 0.05)
    res <- res[pathway %in% keep$mainPathways,]
    res <- res[order(padj),]
    fwrite(res, file = paste0("text/fgsea_BP_", x, ".csv"), sep = ",")
})

