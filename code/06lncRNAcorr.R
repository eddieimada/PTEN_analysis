###########################################################################
### Study gene expression upon PTEN status in TCGA
### Eddie Imada & Luigi Marchionni
###########################################################################
### load libraries
library(ggplot2)
library(reshape2)
library(ggpubr)
rm(list=ls())
### Load objs
### THIS OBJECT IS NOT INCLUDED IN THIS REPOSITORY DUE IT'S LARGE SIZE, BUT CAN BE OBTAINED FROM THE RECOUNT PORTAL
load("~/Dropbox (MechPred)/FANTOM6/TCGA/objs/tcgaEset.rda")


### Rename some stuff for nice plotting
exp <- exprs(TCGAeset)
pheno <- pData(TCGAeset)
feat <- fData(TCGAeset)

SampleTypes <- pheno$cgc_sample_sample_type
SampleTypes <- gsub("Additional - New Primary", "Primary Tumor", SampleTypes)
SampleTypes <- gsub("Additional Metastatic", "Metastatic", SampleTypes)
SampleTypes <- gsub("Additional - New Primary", "Primary Tumor", SampleTypes)

TissueTypes <- pheno$gdc_cases.tissue_source_site.project
TissueTypes <- gsub("Uterine Corpus Endometrial Carcinoma", "Uterine Endometrial Carcinoma", TissueTypes )
TissueTypes <- gsub("Kidney renal clear cell carcinoma", "Kidney RCC carcinoma", TissueTypes )
TissueTypes <- gsub("Pheochromocytoma and Paraganglioma" , "Pheochromocytoma/Paraganglioma" , TissueTypes )
TissueTypes <- gsub("Lymphoid Neoplasm Diffuse Large B-cell Lymphoma", "Diffuse Large B-cell Lymphoma", TissueTypes )
TissueTypes <- gsub("Ovarian serous cystadenocarcinoma", "Ovarian serous adenocarcinoma", TissueTypes )
TissueTypes <- gsub("Cervical squamous cell carcinoma and endocervical adenocarcinoma", "CSCC/Endocervical adenocarcinoma", TissueTypes )
TissueTypes <- gsub("Head and Neck squamous cell carcinoma" , "HNSCC" , TissueTypes )
TissueTypes <- gsub("Ovarian serous cystadenocarcinoma", "Ovarian serous adenocarcinoma", TissueTypes )



### Genes to plot
gnsList <- list(c("CATG00000038715", "ENSG00000186115", "ENSG00000171903", "ENSG00000213903", "ENSG00000111144"),
c("CATG00000079217", "ENSG00000183580", "ENSG00000158715", "ENSG00000167751", "ENSG00000227418", "ENSG00000142515", "ENSG00000157214"),
c("CATG00000000330", "ENSG00000171862"),
c("CATG00000117664", "ENSG00000151025"),
c("CATG00000045713", "ENSG00000225258"))

# Plot correlations
lapply(gnsList, function(gns){
mat <- exp[gns,]
nmsMat <- feat[gns,]$geneName[feat[gns,]$geneID == rownames(mat)]
names(nmsMat) <- gns


### Split by tissue
mat.split <- split(data.frame(log2(t(mat)+1)), TissueTypes)


corrList <- lapply(2:nrow(mat), function(i){
    corr <- unlist(lapply(mat.split, function(x){
        cor(x[,1], x[,i])
    }))
    p.value <- unlist(lapply(mat.split, function(x){
        cor.test(x[,1], x[,i])$p.value
    }))
    cbind.data.frame(corr, p.value)
})
names(corrList) <- rownames(mat)[2:nrow(mat)]
matcorr <- do.call("cbind", corrList)
write.csv(matcorr, file=paste0("text/09", rownames(mat)[1],"_corr.csv"))

dflog <- cbind.data.frame(log2(t(mat)+1), TissueTypes, SampleTypes)
dflog <- dflog[!is.na(dflog$SampleTypes), ]
dflog <- dflog[dflog$SampleTypes != "Recurrent Tumor", ]
dflog$SampleTypes <- gsub("Primary Blood Derived Cancer - Peripheral Blood", "Primary Tumor", dflog$SampleTypes)
refgene <- names(dflog)[1]
sapply(names(dflog)[c(-1,-ncol(dflog))], function(name){
    p <- ggplot(dflog, aes_string(x=refgene, y=name))+
        geom_point(size = 0.5, alpha = 0.5, aes(colour = factor(SampleTypes))) +
        scale_colour_manual(values=c("green", "black","deeppink1")) +
        geom_smooth(method=lm, se = FALSE) +
        facet_wrap(~TissueTypes, nrow = 8, scales = "fixed") +
        theme_bw() +
        xlab(paste0(refgene," / ", nmsMat[refgene], " - Expression (log2)")) +
        ylab(paste0(name," / ", nmsMat[name], " - Expression (log2)")) +
        stat_cor(
            method = "pearson", label.sep = "\n",
            label.x.npc = "left", label.y.npc = "top", size = 2.8, color = "red"
        )
    
    ggsave(p, file=paste0("figs/09", refgene, "-",nmsMat[name],"_corrByTissue.jpeg"), height = 15, width = 12)
})
rm(mat, mat.split, matcorr,dflog)
gc()
})
