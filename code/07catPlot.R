### Clean
rm(list=ls())
library(matchBox)
load("~/Dropbox (MechPred)/Projects/PtenERG/Manuscript/objs/DGE_prad.rda")

xde <- read.csv("~/Dropbox (MechPred)/Projects/PtenERG/XDE/text/PTENposVSneg.csv", stringsAsFactors = F)
xde <- xde[,c(1,6)]
names(xde) <- c("symbol", "value")
xde$value <- xde$value * -1
tcga <- cbind.data.frame(symbol=tGnr$PTEN_NEGvsPTEN_POS$geneName, value=tGnr$PTEN_NEGvsPTEN_POS$t)
tcga <- filterRedundant(tcga)
tcga <- tcga[!is.na(tcga$symbol),]
allData <- list(xde,tcga)

mergedDf <- mergeData(allData, idCol = 1, byCol = 2)
HypPI <- calcHypPI(mergedDf, expectedProp = NULL)
CAT_UP <- computeCat(mergedDf, size = 500)
CAT_DN <- computeCat(mergedDf, size=500, decreasing = FALSE)
CAT <- list(Up=CAT_UP$.value.vs..value.1, 
            Down=CAT_DN$.value.vs..value.1)

png("./figs/CAT.png", width=2200, height = 1100, res=330)
plotCat(CAT, preComputedPI = HypPI, main = "Correspondance-At-the-Top between Up- and Down-regulated genes", maxYlim = 0.3, spacePts = 50, col = c("steelblue", "goldenrod1"), pch = 19, lty = 1, lwd = 1.5 )
dev.off()

xde <- read.csv("~/Dropbox (MechPred)/Projects/PtenERG/XDE/text/PTENinERGPosxde.csv", stringsAsFactors = F)
xde <- xde[,c(1,6)]
names(xde) <- c("symbol", "value")
xde$value <- xde$value * -1
tcga <- cbind.data.frame(symbol=tGnr$PTEN_NEGvsPTEN_POSinERGpos$geneName, value=tGnr$PTEN_NEGvsPTEN_POS$t)
tcga <- filterRedundant(tcga)
tcga <- tcga[!is.na(tcga$symbol),]
allData <- list(xde,tcga)

mergedDf <- mergeData(allData, idCol = 1, byCol = 2)
HypPI <- calcHypPI(mergedDf, expectedProp = NULL)
CAT_UP <- computeCat(mergedDf, size = 500)
CAT_DN <- computeCat(mergedDf, size=500, decreasing = FALSE)
CAT <- list(Up=CAT_UP$.value.vs..value.1, 
            Down=CAT_DN$.value.vs..value.1)

png("./figs/CAT.png", width=2200, height = 1100, res=330)
plotCat(CAT, preComputedPI = HypPI, main = "Correspondance-At-the-Top between Up- and Down-regulated genes", maxYlim = 0.3, spacePts = 50, col = c("steelblue", "goldenrod1"), pch = 19, lty = 1, lwd = 1.5 )
dev.off()