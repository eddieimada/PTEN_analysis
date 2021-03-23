### Loading libraries
library(biomaRt)
library(XDE)
library(GeneMeta)
library(siggenes)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
### Cleaning enviroment and setting wd
rm(list=ls())
setwd("~/Dropbox (MechPred)/Projects/PtenERG/XDE/")

### Getting annotation from bioMart
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ann <- getBM(attributes = c("hgnc_symbol", "description"),
             mart = mart)
ann$description <- gsub(" \\[Source.+", "", ann$description)


### Loading esets
load("~/Dropbox (MechPred)/Projects/PtenERG/TCGA/objs/tcga_mrnaEset.rda")
tcga <- mrnaTCGAeset

load("~/Dropbox (MechPred)/Projects/PtenERG/NaturalHistory/objs/eset.rda")
natHist <- eset

load("~/Dropbox (MechPred)/Projects/PtenERG/hpfsPtenERG/objs/hpfsEset.rda")
hpfs <- hpfsEset

### Getting and keeping only genes in common across all cohorts
TCGAids <- rownames(exprs(tcga))
NatHistids <- rownames(exprs(natHist))
HPFSids <- rownames(exprs(hpfs))

keep <- Reduce(intersect, list(TCGAids,NatHistids,HPFSids))

TCGAmatrix <- exprs(tcga)[keep,]
annTCGA <- featureData(tcga)@data[featureData(tcga)$SYMBOL %in% keep,]
rownames(annTCGA) <- annTCGA$SYMBOL
annTCGA <- annTCGA[rownames(TCGAmatrix),]

### For the natHist cohort, if multiple symbols keep only the first and remove duplicates
natHistMatrix <- exprs(natHist)[keep,]
class(natHistMatrix) <- "numeric"
natHistMatrix <- natHistMatrix[rownames(annTCGA),]

### Creating dicotomous variable for XDE comparison
phenoNatHist <- phenoData(natHist)@data
phenoNatHist$status <- paste(phenoNatHist$PTEN_LOSS, phenoNatHist$ERG, sep = ".")
# phenoNatHist$status <- gsub("\\..+", "", phenoNatHist$status)
# phenoNatHist$status <- gsub("PTEN_POSITIVE", "1", phenoNatHist$status)
# phenoNatHist$status <- gsub("PTEN_NEGATIVE", "0", phenoNatHist$status)
phenoNatHist$status <- gsub("PTEN_NEGATIVE.ERG_NEGATIVE", "1", phenoNatHist$status)
phenoNatHist$status <- gsub("PTEN_NEGATIVE.ERG_POSITIVE", "x", phenoNatHist$status)
phenoNatHist$status <- gsub("PTEN_POSITIVE.ERG_NEGATIVE", "0", phenoNatHist$status)
phenoNatHist$status <- gsub("PTEN_POSITIVE.ERG_POSITIVE", "x", phenoNatHist$status)
phenoNatHist <- phenoNatHist[phenoNatHist$status %in% c(0,1),]
natHistMatrix <- natHistMatrix[,rownames(phenoNatHist)]
phenoNatHist$status <- as.numeric(phenoNatHist$status)
esetNatHist <- ExpressionSet(assayData = natHistMatrix,
                             phenoData = AnnotatedDataFrame(phenoNatHist),
                             featureData = AnnotatedDataFrame(annTCGA))

validObject(esetNatHist)

### Keeping genes in common 
HpfsMatrix <- exprs(hpfs)[keep,]
HpfsMatrix <- HpfsMatrix[rownames(annTCGA),]
class(HpfsMatrix) <- "numeric"


### Creating dicotomous variable for XDE comparison
phenohpfs <- phenoData(hpfs)@data
keep1 <- pData(hpfs)$status == "TUMOR"
phenohpfs <- phenohpfs[keep1,]
HpfsMatrix <- HpfsMatrix[,keep1]
phenohpfs$status <- paste(phenohpfs$pten_pos, phenohpfs$erg_pos, sep = ".")
phenohpfs$status <- gsub("POS", "POSITIVE", phenohpfs$status)
phenohpfs$status <- gsub("NEG", "NEGATIVE", phenohpfs$status)
phenohpfs$status <- gsub("PTEN_NEGATIVE.ERG_NEGATIVE", "1", phenohpfs$status)
phenohpfs$status <- gsub("PTEN_NEGATIVE.ERG_POSITIVE", "x", phenohpfs$status)
phenohpfs$status <- gsub("PTEN_POSITIVE.ERG_NEGATIVE", "0", phenohpfs$status)
phenohpfs$status <- gsub("PTEN_POSITIVE.ERG_POSITIVE", "x", phenohpfs$status)
phenohpfs <- phenohpfs[phenohpfs$status %in% c(0,1),]
HpfsMatrix <- HpfsMatrix[,rownames(phenohpfs)]
phenohpfs$status <- as.numeric(phenohpfs$status)

esetHPFS <- ExpressionSet(assayData = HpfsMatrix,
                          phenoData = AnnotatedDataFrame(phenohpfs),
                          featureData = AnnotatedDataFrame(annTCGA))

validObject(esetHPFS)

### Creating XDE object ExpressionSetList
esetList <- new("ExpressionSetList", .Data=list(esetNatHist,esetHPFS))

### Checking if all datasets contain dicotomous variable for comparisons
all(sapply(esetList, function(x, label){ label %in% varLabels(x)}, label="status"))

### Creating XdeParameter class with deltag model (assumes same deltas for each study)
params <- new("XdeParameter", esetList=esetList, phenotypeLabel="status", one.delta = TRUE)

### Fitting bayesian hierarchical model
# Using empirical values for starting the chain
empirical <- empiricalStart(esetList, phenotypeLabel="status")
firstMcmc(params) <- empirical
# Number of bootstraps
iterations(params) <- 1000
# Saving log files
burnin(params) <- FALSE
# Saving only parameters not indexed by gene or study
output(params)[c("potential", "acceptance", "diffExpressed", "nu",
                 ##"DDelta", 
                 ##"delta", 
                 "probDelta", 
                 ##"sigma2",
                 "phi")] <- 0
# Save only every 2 iterations
thin(params) <- 2
# Directory of the log files
directory(params) <- "logFilesPTENbyERGpos"
# Fit the BH model
xmcmc <- xde(params, esetList)

# Save object
save(xmcmc, file = "./objs/mcmcDeltaTRUE_PTENbyERGneg.rda")
#load("./objs/mcmcDeltaFALSE_PTENbyERGpos.rda")

### Extracting and plotting parameters that are not indexed by gene and platform
getLogs <- function(object){
    params <- output(object)[output(object) == 1]
    params <- params[!(names(params) %in% c("nu", "phi", "DDelta", "delta", "sigma2", "diffExpressed"))]
    names(params)
}
param.names <- getLogs(xmcmc)
params <- lapply(lapply(as.list(param.names), 
                        function(name, object) eval(substitute(object$NAME_ARG, 
                                                               list(NAME_ARG=name))), object=xmcmc), as.ts)
names(params) <- param.names
tracefxn <- function(x, name) plot(x, plot.type="single", col=1:ncol(x), ylab=name)
mapply(tracefxn, params, name=names(params))


### Calculating Bayesian Effect Size (posterior mean of the standardized offsets)
bayesianEffectSize(xmcmc) <- calculateBayesianEffectSize(xmcmc)
### Calculating the posterior average for indicators of concordant and discordant differential expression
posteriorAvg(xmcmc) <- calculatePosteriorAvg(xmcmc, NDIFF = 2)

# Save object
save(xmcmc, file = "./objs/mcmcDeltaTRUE_PTENbyERGneg.rda")
#load("./objs/mcmcDeltaFALSE_PTENbyERGpos.rda")
#load("./objs/esetListPTENposXPTENneginERGpos.rda")


#Extract BES and PP
BES <- bayesianEffectSize(xmcmc)
postAvg <- posteriorAvg(xmcmc)
postAvg <- postAvg[,c(1,2)]
all(rownames(postAvg) == rownames(BES))
merge <- cbind(postAvg, BES)
merge <- cbind(merge, (merge[,3] + merge[,4]))
colnames(merge)[3:5] <- c("BES_NatHist", "BES_HPFS", "BES_SUM")
sum <- merge[,1] - merge[,2]
sum[sum < 0] <- 0
mergedOrderedbyConc <- merge[order(-merge[,1], -merge[,5] ), ]
mergedOrderedbyES <- merge[order( -merge[,5], -merge[,1] ), ]
postAvgES <- mergedOrderedbyES[,c(1,2)]
BESES <- mergedOrderedbyES[,c(3,4)]

#Plot Figure
png("figs/PTENinERGNegative/postAvgBES_PTENbyERGneg_byES_UP.png", width=1500, height = 2500, res=300)

BESflip <- rbind(head(BESES,25), tail(BESES,25))
colnames(BESflip) <- c("Natural History", "HPFS/HPS")
BES.m <- melt(BESflip)
names(BES.m)[2] <- "Cohort" 
plot2 <- ggplot(BES.m, aes(Var1, value,fill = Cohort)) +
    geom_bar(position="identity", stat = "identity",alpha = 0.5) +
    coord_flip() +
    xlab("") +
    ylab("Bayesian Effect Size") +
    theme_bw() +
  theme(legend.position = "bottom")
  
plot2
dev.off()


write.csv(merge, file = "./text/PTENinERGnegxde.csv", row.names = TRUE)
