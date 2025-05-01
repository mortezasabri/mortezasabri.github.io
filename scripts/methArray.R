# Morteza Sabri
# 13.01.2025


# required packages
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
pkgs <- c("ggplot2", "dplyr", "RColorBrewer", "limma", "minfi", 
          "methylationArrayAnalysis",
          "IlluminaHumanMethylationEPICv2manifest",
          "IlluminaHumanMethylationEPICv2anno.20a1.hg38", 
          "TxDb.Hsapiens.UCSC.hg38.knownGene", "DMRcatedata", 
          "EnhancedVolcano", "karyoploteR", "BSgenome.Hsapiens.UCSC.hg38", 
          "BSgenome.Hsapiens.UCSC.hg38.masked") # humarray
new.pkg <- pkgs[!(pkgs %in% utils::installed.packages()[, "Package"])]
if (length(new.pkg)) BiocManager::install(new.pkg, update = FALSE, force = TRUE)
if (!requireNamespace("methyPre")) devtools::install_github("metamaden/methyPre")


write.txt.xlsx <- function(data, outdir, prefix = deparse(substitute(data)) ) {
    utils::write.table(as.data.frame(data), 
                       file = paste0(outdir, prefix, ".txt"), 
                       quote = F, sep = "\t", row.names = FALSE)
    writexl::write_xlsx(data.frame(data), paste0(outdir, prefix,".xlsx"))
}

# ------- Variables
project <- "ppdT1_vs_ppdT0"
normMethod <- "preprocessQuantile" # preprocessQuantile or preprocessFunnorm
gnom <- "hg38"     # c("hg19", "hg38", "mm10")
dataType <- "array"     # c("array", "sequencing")
arrayType <- "EPICv2"     # c("EPICv2", "EPICv1", "EPIC", "450K")
betaCutOff <- 0.1
pCutoff <- 0.05
palettes <- pal <- RColorBrewer::brewer.pal(8,"Dark2")

# sub_percent_data <- 10  # 1/10 subset the data (rows)

partent_outdir <- out <- paste0("out/", project, "/")
if(!dir.exists(partent_outdir)) dir.create(partent_outdir, recursive = TRUE)

idat <- paste0( "idat/", project, "/")
if(!dir.exists(idat)) dir.create(idat, recursive = TRUE)

input <- "dat"

grps <- c("condition", "Time")  # groups (factors) based of the names of gr object

var1 <- "condition"
condition1 <- "Healthy"
condition1_abrv <- "H"
condition2 <- "PPD"
condition2_abrv <- "D"

var2 <- "Time"
time1 <- "pregnant"
time1_abrv <- "T0"
time2 <- "post"
time2_abrv <- "T1"


if (project == "ppdT1_vs_ppdT0") {
        sbs <- "PPD"
        vars <- "condition"
        grps1 <- "pregnant"
        grps2 <- "post"
        projectTitle <- "PPD_week_12 vs. PPD_day_2"
    } else if (project == "healthyT1_vs_healthyT0") {
        sbs <- "Healthy"
        vars <- "condition"
        grps1 <- "pregnant"
        grps2 <- "post"
        projectTitle <- "healthy_week_12 vs. healthy_day_2"
    } else if (project == "ppdT0_vs_healthyT0") {   #new
        sbs <- "pregnant"
        vars <- "Time"
        grps1 <- "Healthy"
        grps2 <- "PPD"
        projectTitle <- "PPD_day_2 vs. healthy_day_2"
    } else if (project == "ppdT1_vs_healthyT1") {   #new
        sbs <- "post"
        vars <- "Time"
        grps1 <- "Healthy"
        grps2 <- "PPD"
        projectTitle <- "PPD_week_12 vs. healthy_week_12"
    } else if (project == "t1_vs_t0") {
        if (exists("sbs")) rm(sbs)
        if (exists("vars")) rm(vars)
        grps1 <- "pregnant"
        grps2 <- "post"
        fillOfChoice <- "Time"
        projectTitle <- "week_12 vs. day_2"
    } else if (project == "ppd_vs_healthy") {
        if (exists("sbs")) rm(sbs)
        if (exists("vars")) rm(vars)
        grps1 <- "Healthy"
        grps2 <- "PPD"
        fillOfChoice <- "condition"
        projectTitle <- "PPD vs. healthy"
    } else if (project == "ppdT1_vs_others") {
        if (exists("sbs")) rm(sbs)
        if (exists("vars")) rm(vars)
        grps1 <- "pregnant"
        grps2 <- "post"
        fillOfChoice <- "Time"
        projectTitle <- "PPD_week_12 vs. others"
}

# fill object for ggplot2, the one that it's gonna be in the plots 
if (exists("vars")) {
    ifelse(vars == "condition", fill <- "Time", fill <- "condition")
} else {
    fill <- grps
}

if (!exists("fillOfChoice")) {
    fillOfChoice <- fill[1]   # in case of just one input
}


# minifi pipeline ##############################################################

# targets <- minfi::read.metharray.sheet(input, pattern="dat/SampleSheet.csv")
targets0 <- read.csv("dat/targets.csv", row.names = 1)
str(targets0)
targets0 <- as.data.frame(lapply(targets0, as.factor))
str(targets0)
targets0$condition <- relevel(targets0$condition, "Healthy")
targets0$Time <- relevel(targets0$Time, "pregnant")


if (exists("sbs")) {
    ind <- which(colnames(targets0) == vars)
    targets <- targets0[targets0[, ind] == sbs, ]
} else {
    targets <- targets0
}


## ----- reading data ----
rgSet <- minfi::read.metharray.exp(targets = targets) #base = input)
if (exists("sub_percent_data")) {
    keep_probes <- sample(rownames(rgSet), nrow(rgSet)/10)
    rgSet <- rgSet[keep_probes, ]
}

minfi::annotation(rgSet)
# Check sample pairing
minfi::sampleNames(rgSet)



v <- paste(targets$Sample_Name, targets[, grps[1]], targets[, grps[2]], sep = "_")
v <- sub(condition1, condition1_abrv, v)
v <- sub(time1, time1_abrv, v)
v <- sub(condition2, condition2_abrv, v)
v <- sub(time2, time2_abrv, v)
targets$ID <- v

# gr object
gr <- data.frame(targets[, c(grps[1], grps[2]), drop = FALSE], row.names = targets$ID)


minfi::sampleNames(rgSet) <- targets$ID
Biobase::pData(rgSet)$id <- targets$ID
Biobase::pData(rgSet)$group <- gr
rgSet   # weird that there is no colnames



# Quality control --------------------------------------------------------------
outdir <- paste0(partent_outdir, "1_qc/")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# calculate the detection p-values
detP <- minfi::detectionP(rgSet)


# examine mean detection p-values across all samples to identify any failed samples
df <- data.frame(value = colMeans(detP), gr)
g <- ggplot2::ggplot(df, ggplot2::aes(x = factor(rownames(df), 
                                                 levels = rownames(df)), 
                                      y = value, 
                                      fill = df[, fillOfChoice] )) +
     ggplot2::geom_bar(stat='identity') +
     ggplot2::labs(x = "sample", 
                   y = "p value", 
                   title = "Mean detection p-values across all samples") +
     ggplot2::theme_classic() + 
     ggplot2::theme(legend.position = "none", 
                    axis.text.x = ggplot2::element_text(angle = 90, 
                                                        vjust = 0.5, 
                                                        hjust=1))
ggplot2::ggsave(paste0(outdir, "qc_avgPval_barPlot.png"), 
                 g, dpi = 200, width = 10, height = 6)

if (!exists("sub_percent_data")) {
    minfi::qcReport(rgSet, sampNames=targets$ID, sampGroups = targets[, fill[1]], 
                pdf = paste0(outdir, "/qc_rgSet_qcReport.pdf"))
}



# remove poor quality samples
mean(keep <- colMeans(detP) < 0.05)
rgSet <- rgSet[,keep]
table(keep)
targets <- targets[keep,] # remove poor quality samples from targets data
detP <- detP[,keep] # remove poor quality samples from detection p-value table
dim(detP)




# Normalization ----------------------------------------------------------------
outdir <- paste0(partent_outdir, "2_norm/")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# normalize the data
mSetSq <- tryCatch({
    if (normMethod == "preprocessQuantile") {
        mSetSq <- minfi::preprocessQuantile(rgSet)
    } else if (normMethod == "preprocessFunnorm") {
        mSetSq <- minfi::preprocessFunnorm(rgSet)
    } else {
        stop("define the normMethod obj")
    }
}, error = function(e) {
    # Code to run if an error occurs
    print("An error occurred! Using preprocessNoob() ...")
    preprocessNoob(rgSet)
})


# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)


# visualize what the data looks like before and after normalization
if (length(fill) >= 1) {
    grDevices::png(paste0(outdir, "norm_", "beforeAfterNorm_linePlot.png"), 
                   units = "in", height = 8, width = 12, res = 300)
    par(mfrow=c(1, 2))
    densityPlot(rgSet, sampGroups = targets[, fill[1]], main="Raw", legend=FALSE)
    graphics::legend("top", legend = levels(factor(targets[, fill[1]])), 
           text.col= palettes)
    densityPlot(getBeta(mSetSq), sampGroups = targets[, fill[1]],
                main="Normalized", legend=FALSE)
    graphics::legend("top", legend = levels(factor(targets[, fill[1]])), 
           text.col= palettes)
    dev.off()
} 

if (length(fill) >= 2) {
    grDevices::png(paste0(outdir, "norm_", "beforeAfterNorm_linePlot2.png"), 
                   units = "in", height = 8, width = 12, res = 300)
    par(mfrow=c(1, 2))
    densityPlot(rgSet, sampGroups = targets[, fill[2]], main="Raw", legend=FALSE)
    graphics::legend("top", legend = levels(factor(targets[, fill[2]])), 
           text.col = palettes)
    densityPlot(getBeta(mSetSq), sampGroups = targets[, fill[2]],
                main="Normalized", legend=FALSE)
    graphics::legend("top", legend = levels(factor(targets[, fill[2]])), 
           text.col = palettes)
    dev.off()
} 



# MDS
n <- pData(mSetRaw)$id
g <- pData(mSetRaw)$group
if (length(fill) >= 1 ){
    grDevices::png(paste0(outdir, "norm_", "rawMethylSet_mdsPlot.png"), 
                   units = "in", height = 8, width = 8, res = 300)
    minfi::mdsPlot(mSetRaw, sampNames = n, sampGroups = g[, fill[1]], 
                   main = "MDS")
    dev.off()
} else if (length(fill) == 2){
    grDevices::png(paste0(outdir, "norm_", "rawMethylSet_mdsPlot2.png"), 
                   units = "in", height = 8, width = 8, res = 300)
    minfi::mdsPlot(mSetRaw, sampNames = n, sampGroups = g[, fill[2]], 
                   main = "MDS")
    dev.off()
} else {
    message("minfi::mdsPlot printing just first two fill object")
}




# Filtering --------------------------------------------------------------------
outdir <- paste0(partent_outdir, "3_filt/")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq), rownames(detP)), ] 

# remove any probes that have failed in one or more samples
dim(mSetSq)
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
mSetSqFlt <- mSetSq[keep,]
dim(mSetSqFlt)




# if your data includes males and females, remove probes on the sex chromosomes
anno <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::IlluminaHumanMethylationEPICv2anno.20a1.hg38
anno <- minfi::getAnnotation(anno, what = "everything", lociNames = NULL,
                             orderByLocation = FALSE, dropNonMapping = FALSE)

keep <- !(featureNames(mSetSqFlt) %in% anno$Name[anno$chr %in% 
                                                        c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
dim(mSetSqFlt)

# remove probes with SNPs at CpG site
mSetSqFlt <- minfi::dropLociWithSnps(mSetSqFlt)
dim(mSetSqFlt)




# exclude cross reactive probes 
sum(duplicated(sub("_.*$", "", featureNames(mSetSqFlt))))
keep <- !(sub("_.*$", "", featureNames(mSetSqFlt)) %in% methyPre::pidsley.crxcg)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,] 
dim(mSetSqFlt)



# calculate M-values for statistical analysis
# Assuming `myData` is an RGChannelSet, MethylSet, or GenomicRatioSet object
mVals <- minfi::getM(mSetSqFlt)
head(mVals[,1:5])
# saveRDS(mVals, paste0(idat, "mVals.rds") )

bVals <- minfi::getBeta(mSetSqFlt)
head(bVals[,1:5])
# saveRDS(bVals, paste0(idat, "bVals.rds") )


dBetaDF <- data.frame(cpg = rownames(bVals) ,
                      gr2 = rowMeans(bVals[, gr[, fillOfChoice] == grps2]), 
                      gr1 = rowMeans(bVals[, gr[, fillOfChoice] == grps1]))
dBeta <- dBetaDF$gr2 - dBetaDF$gr1
dBeta <- data.frame(cpg = dBetaDF$cpg, beta = dBeta, 
                        row.names = dBetaDF$cpg)
dBetaCut <- subset(dBeta, beta > betaCutOff)
dim(dBetaCut)



# PCA
# Perform PCA on the transposed beta values (samples as rows, probes as columns)
pca_result <- prcomp(t(bVals), scale. = TRUE)
# Check PCA results
summary(pca_result)  # Variance explained by each principal component
# Extract PCA scores for the first two principal components
df <- as.data.frame(pca_result$x[, 1:2])
colnames(df) <- c("PC1", "PC2")
# Add sample metadata (e.g., group information) to the PCA scores
# Replace 'Group' with your metadata column
df <- cbind(df, pData(mSetSqFlt)[, fill, drop = FALSE])

# Plot the PCA with ggplot2
if (length(fill) == 1) {
    if(exists("g")) rm(g)
    g <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, 
                                                         color = df[, fill[1]], 
                                                  label = rownames(df))) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
        ggplot2::labs(title = "PCA of Methylation Data: Filtered Beta values",
                      x = paste0("PC1 (", 
                                 round(summary(pca_result)$importance[2, 1] * 100,
                                       1), "% Variance)"),
                      y = paste0("PC2 (", 
                                 round(summary(pca_result)$importance[2, 2] * 100,
                                                1), "% Variance)")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.title = ggplot2::element_blank())
} else if (length(fill) >= 2) {
    if(exists("g")) rm(g)
    
    if (which(fill == fillOfChoice) == 1) {
        g <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, 
                                                             color = df[, fill[1]], 
                                                             shape = df[, fill[2]],
                                                      label = rownames(df))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
            ggplot2::labs(title = "PCA of Methylation Data: Filtered Beta values",
                          x = paste0("PC1 (", 
                                     round(summary(pca_result)$importance[2, 1] * 100,
                                           1), "% Variance)"),
                          y = paste0("PC2 (", 
                                     round(summary(pca_result)$importance[2, 2] * 100,
                                           1), "% Variance)")) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.title = ggplot2::element_blank())
    } else if (which(fill == fillOfChoice) == 2) {
        g <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, 
                                              color = df[, fill[2]], 
                                              shape = df[, fill[1]],
                                              label = rownames(df))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
            ggplot2::labs(title = "PCA of Methylation Data: Filtered Beta values",
                          x = paste0("PC1 (", 
                                     round(summary(pca_result)$importance[2, 1] * 100,
                                           1), "% Variance)"),
                          y = paste0("PC2 (", 
                                     round(summary(pca_result)$importance[2, 2] * 100,
                                           1), "% Variance)")) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.title = ggplot2::element_blank())
    }
} else {
    rm(g)
    warning("No fill object (no grouping)")
}
ggplot2::ggsave(paste0(outdir, "qc_bValsFilt_pcaPlot.png"), 
                g, dpi = 200, width = 10, height = 6)





# MDS
if (length(fill) >= 1) {
    grDevices::png(paste0(outdir, "filt_", "mValsFilt_mdsPlot.png"), 
                   units = "in", height = 8, width = 10, res = 300)
    limma::plotMDS(mVals, top=1000, gene.selection="common", 
                   col=pal[factor(targets[, fill[1]])], cex=0.5)
    graphics::legend("right", legend=levels(factor(targets[, fill[1]])), text.col=pal,
           cex=0.65, bg="white")
    dev.off()
} 
if (length(fill) == 2) {
    grDevices::png(paste0(outdir, "filt_", "mValsFilt_mdsPlot2.png"), 
                   units = "in", height = 8, width = 10, res = 300)
    limma::plotMDS(mVals, top = 1000, gene.selection="common", 
                   col= pal[factor(targets[, fill[2]])], cex=0.5)
    graphics::legend("right", legend=levels(factor(targets[, fill[2]])), text.col=pal,
           cex=0.65, bg= "white")
    dev.off()
} else {
    message("limma::plotMDS printing just first two fill object")
}




# density plot
if (length(fill) >= 1) {
    grDevices::png(paste0(outdir, "filt_", "mVals_bValsFilt_densityPlot.png"), 
                   units = "in", height = 8, width = 14, res = 300)
    par(mfrow = c(1,2))
    minfi::densityPlot(bVals, sampGroups=targets[, fill[1]], main="Beta values", 
            legend=FALSE, xlab="Beta values")
    graphics::legend("top", legend = levels(factor(targets[, fill[1]])), 
       text.col= palettes)
    minfi::densityPlot(mVals, sampGroups=targets[, fill[1]], main="M-values", 
            legend=FALSE, xlab="M values")
    graphics::legend("topleft", legend = levels(factor(targets[, fill[1]])), 
            text.col = palettes)
    dev.off()
} 
if (length(fill) == 2) {
    grDevices::png(paste0(outdir, "filt_", "mVals_bValsFilt_densityPlot2.png"), 
                   units = "in", height = 8, width = 14, res = 300)
    par(mfrow = c(1,2))
    minfi::densityPlot(bVals, sampGroups=targets[, fill[2]], main="Beta values", 
                       legend=FALSE, xlab="Beta values")
    graphics::legend("top", legend = levels(factor(targets[, fill[2]])), 
                     text.col= palettes)
    minfi::densityPlot(mVals, sampGroups=targets[, fill[2]], main="M-values", 
                       legend=FALSE, xlab="M values")
    graphics::legend("topleft", legend = levels(factor(targets[, fill[2]])), 
                     text.col = palettes)
    dev.off()
} else if (length(fill) > 2) {
    message("minfi::densityPlot printing just first two fill object")
}   

# save.image(paste0(idat, "filtData.RData"))

# DE -------------------------
outdir <- paste0(partent_outdir, "4_de/")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# design matrix


# use the above to create a design matrix
# design <- model.matrix(~ Subject + Time*condition, data = targets)
# design <- model.matrix(~ 0 + Subject + Time*condition, data = targets)

if (length(fill) == 1) {
    gr_oi1 <- targets[, fill[1]]
} else if (length(fill) == 2) {
    gr_oi1 <- targets[, fill[1]]
    gr_oi2 <- targets[, fill[2]]
} else if (length(fill) == 3) {
    gr_oi1 <- targets[, fill[1]]
    gr_oi2 <- targets[, fill[2]]
    gr_oi3 <- targets[, fill[3]]
} else {
    message("length of fill object is more than 3")
}




# ----------- design ------------
if (length(levels(gr_oi1)) == 2) {
    design <- model.matrix(~ gr_oi1, data = targets)
    n <- paste0(levels(gr_oi1)[2], "_vs_", levels(gr_oi1)[1])
    colnames(design)[ncol(design)] <- n
    names(attributes(design)$contrasts) <- n
    
    coef <- ncol(design)
} else {
    stop("re-write the codes from now on based on 3 or more levels")
}

if (exists("gr_oi2")) {
    if (project == "t1_vs_t0") {
        design <- model.matrix(~ gr_oi2, data = targets)
        n <- paste0(levels(gr_oi2)[2], "_vs_", levels(gr_oi2)[1])
        colnames(design)[ncol(design)] <- n
        names(attributes(design)$contrasts) <- n
        
        coef <- ncol(design)
    } else if (project == "ppd_vs_healthy") {
        design <- model.matrix(~ gr_oi1, data = targets)
        n <- paste0(levels(gr_oi1)[2], "_vs_", levels(gr_oi1)[1])
        colnames(design)[ncol(design)] <- n
        names(attributes(design)$contrasts) <- n
        
        coef <- ncol(design)
    } else if (project == "ppdT1_vs_others") {
        design <- model.matrix(~ gr_oi1 + gr_oi2 + gr_oi1*gr_oi2, data = targets)
        n <- paste0(levels(gr_oi1)[2], "_vs_", levels(gr_oi1)[1])
        colnames(design)[2] <- n
        names(attributes(design)$contrasts)[1] <- n
        
        n <- paste0(levels(gr_oi2)[2], "_vs_", levels(gr_oi2)[1])
        colnames(design)[3] <- n
        names(attributes(design)$contrasts)[2] <- n
        
        colnames(design)[4] <- paste0(colnames(design)[2],":", colnames(design)[3])
        
        coef <- ncol(design)
    } else if (length(levels(gr_oi1)) == 2 & length(levels(gr_oi2)) == 2) {
        design <- model.matrix(~ gr_oi1 + gr_oi2 + gr_oi1*gr_oi2, data = targets)
        n <- paste0(levels(gr_oi1)[2], "_vs_", levels(gr_oi1)[1])
        colnames(design)[ncol(design)] <- n
        names(attributes(design)$contrasts) <- n
        
        
        coef <- ncol(design)
    } else {
        stop("re-write the codes from now on based on 3 or more levels")
    }
}

fit <- limma::lmFit(mVals, design)
fit2 <- limma::eBayes(fit)

# look at the numbers of DM CpGs at FDR < 0.05
summary(limma::decideTests(fit2, p.value = pCutoff))
grDevices::png(paste0(outdir, "sigMultiTest_summaryPlot.png"), units = "in", 
               height = 8, width = 8, res = 300)
plot(summary(limma::decideTests(fit2)))
dev.off()

# get the table of results for the first contrast
annoSub <- anno[match(rownames(mVals), anno$Name),
                c(1:4,12:19,24:ncol(anno))]


DMPs <- limma::topTable(fit2, num=Inf, coef = coef, genelist = annoSub, p.value = 1)
if (nrow(DMPs) == 0) stop("No differentially methylated CpGs to continue with")
head(DMPs)
saveRDS(DMPs, paste0(idat, "DMPs.rds") )
write.txt.xlsx(DMPs, outdir)





# region chromosome plot
cpgs <- GRanges(seqnames = DMPs[, "chr"], 
        ranges = IRanges(start = DMPs[, "pos"], 
                         end = DMPs[, "pos"], 
                         names = rownames(DMPs)),
        strand = Rle(strand(DMPs[, "strand"])), 
        pval = DMPs$adj.P.Val, 
        beta2 = dBetaDF[match(rownames(DMPs), rownames(dBetaDF)), "gr2"],
        beta1 = dBetaDF[match(rownames(DMPs), rownames(dBetaDF)), "gr1"])
cpgs$dBeta <- cpgs$beta2 - cpgs$beta1
cpgs$DM_status <- ifelse(cpgs$dBeta > 0, "hyper", "hypo")
cpgs$DM_status <- ifelse(cpgs$pval < pCutoff, cpgs$DM_status, "none")
genome(cpgs) <- gnom
cpgs <- trim(cpgs)  # This ensures all regions fit within chromosome lengths
seqlengths(cpgs) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(cpgs)]
cpgs <- regioneR::toGRanges(cpgs, genome = gnom)

library(karyoploteR)
kp <- karyoploteR::plotKaryotype(genome = gnom)
# karyoploteR::kpPlotRegions(kp, data=cpgs)
kpPlotRegions(kp, data=cpgs, r0=0, r1=0.5)
kpPlotDensity(kp, data=cpgs, r0=0.5, r1=1)




############################################################ vis barplot regions islands
annots <- grep(gnom, annotatr::builtin_annotations(), value = T)
annotations <- annotatr::build_annotations(genome = gnom, annots)
dm_annotated <- annotatr::annotate_regions(cpgs, annotations = annotations)
dm_annsum <- annotatr::summarize_annotations(
    annotated_regions = dm_annotated,
    quiet = TRUE)
print(dm_annsum)


dm_random_regions <- annotatr::randomize_regions(
    regions = cpgs,
    allow.overlaps = TRUE,
    per.chromosome = TRUE)
dm_random_annotated <- annotatr::annotate_regions(
    regions = dm_random_regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = TRUE)

annotatr::plot_annotation(annotated_regions = dm_annotated,
                          annotated_random = dm_random_annotated,
                          annotation_order = annots,
                          plot_title = '# of Sites Tested for DM annotated',
                          x_label = 'knownGene Annotations',
                          y_label = 'Count')
annotatr::plot_coannotations(annotated_regions = dm_annotated,
    annotation_order = annots,
    axes_label = 'Annotations',
    plot_title = 'Regions in Pairs of Annotations')




# scatter plot
df <- subset(dm_annotated, pval < pCutoff)
g <- annotatr::plot_numerical(annotated_regions = df,
    x = 'beta1',
    y = 'beta2',
    facet = 'annot.type',
    facet_order = annots,
    plot_title = 'Region Methylation: Group 2 vs Group 1',
    x_label = 'Group 1',
    y_label = 'Group 2')
g + ggplot2::geom_point(alpha=0.01)



# --------------- barplot islands
x_order <- c('hyper', 'hypo')
fill_order <- c('hg38_cpg_islands',
    'hg38_cpg_shores',
    'hg38_cpg_shelves',
    'hg38_cpg_inter')
# Make a barplot of the data class where each bar
# is composed of the counts of CpG annotations.
g <- annotatr::plot_categorical(annotated_regions = dm_annotated, 
    x='DM_status', fill='annot.type',
    x_order = x_order, fill_order = fill_order, position='stack',
    plot_title = 'DM Status by CpG Annotation Counts',
    legend_title = 'Annotations',
    x_label = 'DM status',
    y_label = 'Count')

g <- annotatr::plot_categorical(annotated_regions = dm_annotated, 
    x='DM_status', fill='annot.type',
    x_order = x_order, fill_order = fill_order, position='fill',
    plot_title = 'DM Status by CpG Annotation Proportions',
    legend_title = 'Annotations',
    x_label = 'DM status',
    y_label = 'Proportion')

g <- annotatr::plot_categorical(
    annotated_regions = dm_annotated, annotated_random = dm_random_annotated,
    x='DM_status', fill='annot.type',
    x_order = x_order, fill_order = fill_order, position='fill',
    plot_title = 'DM Status by CpG Annotation Proportions',
    legend_title = 'Annotations',
    x_label = 'DM status',
    y_label = 'Proportion')


x_order = c('hg38_custom_ezh2',
    'hg38_genes_1to5kb',
    'hg38_genes_promoters',
    'hg38_genes_5UTRs',
    'hg38_genes_exons',
    'hg38_genes_introns',
    'hg38_genes_3UTRs',
    'hg38_genes_intergenic')
# The orders for the fill labels.
fill_order = c('hyper',
    'hypo',
    'none')
g <- annotatr::plot_categorical(
    annotated_regions = dm_annotated, x='annot.type', fill='DM_status',
    x_order = x_order, fill_order = fill_order, position='fill',
    legend_title = 'DM Status',
    x_label = 'knownGene Annotations',
    y_label = 'Proportion')


# Boxplot









# plot the top 4 most significantly differentially methylated CpGs 
if (length(fill) == 1) {
    par(mfrow=c(2,2))
    sapply(rownames(DMPs)[1:4], function(cpg){
        minfi::plotCpg(bVals, cpg=cpg, pheno= gr_oi1, ylab = "Beta values")
    })
} else if (length(fill) == 2) {
    par(mfrow=c(2,2))
    sapply(rownames(DMPs)[1:4], function(cpg){
        minfi::plotCpg(bVals, cpg=cpg, pheno= gr_oi2, ylab = "Beta values")
    })
}


# histogram
grDevices::png(paste0(outdir, "disPval_histogramBase_pval.png"), units = "in", 
               height = 8, width = 8, res = 300)
graphics::hist(DMPs$P.Value, main = "Histogram of pvalues", xlab = "pvalue")
grDevices::dev.off()

g <- ggplot2::ggplot(DMPs, ggplot2::aes(x=P.Value)) + 
    ggplot2::geom_histogram(bins = 40) + 
    ggplot2::labs(title = "P-value histogram")
ggplot2::ggsave(paste0(outdir, "disPval_histogram_pval.png"), g, dpi = 300)


df <- data.frame(cpg = rownames(DMPs), DMPs)
ind <- intersect(dBetaCut$cpg, df$cpg)
DMPsdBeta <- merge(dBetaCut, df, by = "cpg")
rm(df)
saveRDS(DMPsdBeta, paste0(idat, "DMPsdBeta.rds") )
write.txt.xlsx(DMPsdBeta, outdir)




# Volcano plot
# pvalues from DMPs and delta values from dBeta
if ( all(rownames(dBeta) == rownames(DMPs)) ) {
    df <- data.frame(cpg = dBeta$cpg, 
                     detaBeta = dBeta$beta, 
                     pval = DMPs$P.Value, row.names = dBeta$cpg)
    caption <- paste0("total = ", nrow(df), " variables",
                      "    pos delta: ", 
                      sum(df$detaBeta > betaCutOff & df$pval < pCutoff), 
                      "    neg delta: ", 
                      sum(df$detaBeta < -betaCutOff & df$pval < pCutoff), 
                      "    Cutoffs: beta |", betaCutOff, "| & pval ", pCutoff)
    g <- EnhancedVolcano::EnhancedVolcano(
        df, lab = df$cpg,
        title = paste0("Volcano Plot ", projectTitle), 
        x = "detaBeta",
        y = "pval",
        pCutoff = pCutoff,
        FCcutoff = betaCutOff,
        selectLab = DMPsdBeta$cpg,
        boxedLabels = TRUE,
        drawConnectors = T,
        widthConnectors = 0.3,
        colConnectors = "black",
        pointSize = 1,
        colAlpha = 1,
        caption = caption,
        col = c("grey30", 'grey30', 'grey30', "red2")) +
        ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(paste0(outdir, "volcanoPlot.png"), 
                    g, dpi = 300, height = 17, width = 17)
} else {
    stop("merge dBeta and DMPs dataframes")
}




# Nearest genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- GenomicFeatures::genes(txdb)
df <- DMPsdBeta
query_position <- GRanges(seqnames = df$chr, 
                          ranges = IRanges(start = df$pos, 
                                           end = df$pos),
                          strand = Rle(strand(df$strand)))

# Find the closest gene
nearest_gene_index <- nearest(query_position, genes)
nearest_gene <- genes[nearest_gene_index]

# Get Gene ID and Symbol
gene_id <- as.character(nearest_gene$gene_id)
gene_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                     keys = gene_id, column = "SYMBOL", 
                                     keytype = "ENTREZID", 
                                     multiVals = "first")

# Output the closest gene
df2 <- data.frame(nearest_entrezid = gene_id, 
                  nearest_symbol = gene_symbol, 
                  nearest_chr = seqnames(nearest_gene),
                  nearest_start = start(nearest_gene), 
                  nearest_end = end(nearest_gene))

if (nrow(DMPsdBeta) == nrow(df2)) {
    df <- cbind(DMPsdBeta, df2)
    
    # and put the distance in the table for them
    for (i in 1:nrow(df2)) {
        df$distCpGToNearestGene[i] <- df$nearest_start[i] - df$pos[i]
    }
    DMPsdBeta <- df
    saveRDS(DMPsdBeta, paste0(idat, "DMPsdBeta.rds") )
    write.txt.xlsx(DMPsdBeta, outdir)
}





# Differential methylation analysis of regions
myAnnotation <- DMRcate::cpg.annotate(object = mVals, datatype = dataType,
                                      what = "M",
                                      analysis.type = "differential",
                                      design = design,
                                      contrasts = F, arraytype = arrayType,
                                      coef = coef, 
                                      fdr = ifelse(project == "ppd_vs_healthy", 0.2, pCutoff))
saveRDS(myAnnotation, paste0(idat, "myAnnotation.rds") )


DMRs <- DMRcate::dmrcate(myAnnotation, lambda=1000, C=2, betacutoff = betaCutOff)
# saveRDS(DMRs, paste0(idat, "DMRs.rds") )

results.ranges <- DMRcate::extractRanges(DMRs, genome = gnom)
# saveRDS(results.ranges, paste0(idat, "results.ranges.rds") )
# write.txt.xlsx(results.ranges, outdir)






# with dBeta values -------------------------------------------------------------------- in regions is better
kp <- karyoploteR::plotKaryotype(genome = gnom, plot.type=2)






# --------------------------- Visualization ------------------------------------

# set up the grouping variables and colours
groups <- pal[1:length(unique(targets[, fillOfChoice]))]
names(groups) <- levels(factor(targets[, fillOfChoice]))
cols <- groups[as.character(factor(targets[, fillOfChoice]))]
# draw the plot for the top DMR



get.PlotTracks <- function(results.ranges, bVals, cols, 
                           arrayType, gnom, outdir) {
    grDevices::png(paste0(outdir, "dmr_", dmrIndex, "_heatmap.png"), units = "in", 
                   height = 8, width = 14, res = 300)
    DMRcate::DMR.plot(ranges = results.ranges, dmr = dmrIndex, CpGs = bVals, 
                      phen.col = cols, 
                      what = "Beta", arraytype = arrayType, genome = gnom)
    dev.off()
    message(paste0("Region ", dmrIndex, " is done"))
}


for (dmrIndex in 1:nrow(data.frame(results.ranges))) {
    if (dmrIndex > 20) break("dmrIndex exceeds the limit of 50")
    invisible(gc())
    message(paste0("DMRcate::DMR.plot region ", 
                   dmrIndex, " is in progress"))
    tryCatch({
        R.utils::withTimeout(get.PlotTracks(results.ranges, bVals, cols, 
                                            arrayType, gnom, outdir), 
                             timeout = 120)
        
    }, TimeoutException = function(ex) {
        print("Timeout of 60 seconds reached. Function will not be completely run.")
    }, error = function(e) {
        warning(paste0("Region Index ", dmrIndex, " has been skipped (due to NAs)"))
        dev.off()
    } )
    
}

# system("find out/healthyT1_vs_healthyT0/4_de/ -type f -size -1024 -delete")




# ------------------------ customized visualization

customizedPlotTracks <- function(results.ranges, dmrIndex, annoSub, bVals, 
                                 targets, fill, gnom, chrom, pal, outdir) {
    # extract chromosome number and location from DMR results 
    chrom <- as.character(seqnames(results.ranges[dmrIndex]))
    start <- as.numeric(start(results.ranges[dmrIndex]))
    end <- as.numeric(end(results.ranges[dmrIndex]))
    # add 25% extra space to plot
    minbase <- start - (0.25*(end-start))
    maxbase <- end + (0.25*(end-start))
    
    
    annoSubSort <- annoSub[order(annoSub$chr, annoSub$pos),]
    head(annoSubSort)
    
    bValsSort <- bVals[match(annoSubSort$Name, rownames(bVals)),]
    head(bValsSort)
    
    # create genomic ranges object from methylation data
    cpgData <- GRanges(seqnames=Rle(annoSubSort$chr),
                       ranges=IRanges(start=annoSubSort$pos, 
                                      end=annoSubSort$pos),
                       strand=Rle(rep("*", nrow(annoSubSort))),
                       betas=bValsSort)
    # extract data on CpGs in DMR
    cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])
    
    # methylation data track
    methTrack <- Gviz::DataTrack(range=cpgData, groups= targets[, fill],
                                 genome = gnom, chromosome=chrom, 
                                 ylim=c(-0.05,1.05), col=pal,
                                 type=c("a","p"), name="DNA Meth.\n(beta value)",
                                 background.panel="white", legend=TRUE, cex.title=0.8,
                                 cex.axis=0.8, cex.legend=0.8)
    
    
    # DMR position data track
    dmrTrack <- Gviz::AnnotationTrack(start=start, end=end, genome=gnom, 
                                      name="DMR", 
                                      chromosome=chrom, fill="darkred")
    
    
    
    iTrack <- Gviz::IdeogramTrack(genome = gnom, chromosome = chrom, name="")
    
    gTrack <- Gviz::GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
    
    rTrack <- Gviz::UcscTrack(genome=gnom, chromosome=chrom, track="NCBI RefSeq", 
                              from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                              rstarts="exonStarts", rends="exonEnds", gene="name", 
                              symbol="name2", transcript="name", strand="strand", 
                              fill="darkblue",stacking="squish", name="RefSeq", 
                              showId=TRUE, geneSymbol=TRUE)
    
    
    tracks <- list(iTrack, gTrack, rTrack, methTrack, dmrTrack)
    
    sizes <- c(1, 1, 1, 6, 1) # set up the relative sizes of the tracks
    
    grDevices::png(paste0(outdir, "dmr_", dmrIndex, "_map.png"), units = "in", 
                   height = 8, width = 14, res = 300)
    Gviz::plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
                     add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
    dev.off()
    message(paste0("Region ", dmrIndex, " is done"))
}


for (dmrIndex in 1:nrow(data.frame(results.ranges))) {
    if (dmrIndex > 20) break("dmrIndex exceeds the limit of 50")
    invisible(gc())
    message(paste0("Customized Gviz::plotTracks region ", 
                   dmrIndex, " is in progress"))
    tryCatch({
        R.utils::withTimeout(customizedPlotTracks(results.ranges, dmrIndex, annoSub, bVals, 
                                         targets, fill, gnom, chrom, pal, outdir), timeout = 120)
        
    }, TimeoutException = function(ex) {
        print("Timeout of 60 seconds reached. Function will not be completely run.")
        dev.off()
    }, error = function(e) {
        warning(paste0("Region Index ", dmrIndex, " has been skipped (due to NAs)"))
        dev.off()
    } )
}












# ------------ Volcano
DMPs <- limma::topTable(fit2, num=Inf, coef= coef, genelist= annoSub)
EnhancedVolcano::EnhancedVolcano(DMPs,
                lab = rownames(DMPs),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Differential Methylation Volcano Plot")


# density plot
beta_filtered <- bVals[rowMeans(bVals > 0 & bVals < 1, na.rm = TRUE) > 0.05, ]
for (i in 1:ncol(beta)) {
    tryCatch({
        beta[, i] <- champ.BMIQ(beta[, i], design.v, sampleID = colnames(beta)[i])$beta.corr
    }, error = function(e) {
        message(paste("BMIQ failed for sample:", colnames(beta)[i]))
    })
}
myNorm <- ChAMP::champ.norm(beta = beta_filtered, method = "PBC")  # Alternative: "PBC"

ChAMP::champ.QC(beta = bVals[DMPsdBeta$cpg, ], pheno = gr[, fillOfChoice])


# heatmap
topDMPs <- head(DMPs, 50)
pheatmap::pheatmap(bVals[rownames(topDMPs),], scale = "none", cluster_cols = TRUE)

df <- bVals[DMPsdBeta$cpg, ]
pheatmap::pheatmap(df, scale = "none", cluster_cols = TRUE, 
                   show_rownames = FALSE, 
                   filename = paste0(outdir, "heatmap_bVals.png"),
                   main = paste0("Beta values for top ", nrow(df), " CpGs"))


# --------------- the Most variable probs
fitvar <- missMethyl::varFit(mVals, design = design, coef = c(1, coef))

summary(limma::decideTests(fitvar))
topDV <- missMethyl::topVar(fitvar, coef=2)

par(mfrow=c(2,2))
sapply(rownames(topDV)[1:4], function(cpg){
    plotCpg(bVals, cpg = cpg, pheno= gr[, fillOfChoice], 
            ylab = "Beta values")
})



load(paste0("idat/ppdT1_vs_ppdT0/tmp.RData"))
# ----------------------- GO -------------------------
# Get the significant CpG sites at less than 5% FDR with beta more than 0.1
sigCpGs <- DMPsdBeta$Name[DMPsdBeta$adj.P.Val < 0.01]
length(sigCpGs)

# Get all the CpG sites used in the analysis to form the background
all <- limma::topTable(fit2, num=Inf, coef= coef, genelist= annoSub)$Name
# Total number of CpG sites tested
length(all)


par(mfrow=c(1,1))
gst <- missMethyl::gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE, 
                          array.typ = ifelse(arrayType == "EPICv2","EPIC_V2", arrayType))

missMethyl::topGSA(gst, number=10)







# Cell type composition
pData(rgSet)$Slide <- as.numeric(pData(rgSet)$Slide)
# estimate cell counts
cellCounts <- estimateCellCounts(rgSet, compositeCellType = "Blood", 
                                 referencePlatform = "IlluminaHumanMethylationEPIC")






# 
# -------- ChAMP -----------

library(ChAMP)
# ??IlluminaHumanMethylationEPICv2anno.20a1.hg38
# library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
# Anno <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::IlluminaHumanMethylationEPICv2anno.20a1.hg38
# Anno <- minfi::getAnnotation(Anno, what = "everything", lociNames = NULL, 
#                              orderByLocation = FALSE, dropNonMapping = FALSE)
# Anno <- read.table(gzfile("ref/EPICv2.hg19.manifest.tsv.gz")) 
# Anno <- IlluminaHumanMethylationEPICv2manifest::IlluminaHumanMethylationEPICv2manifest
myLoad <- ChAMP::champ.load("dat2", arraytype = "EPIC")
myLoad <- ChAMP::champ.load("dat2", arraytype = "EPIC", method = "minfi")
CpG.GUI(arraytype="EPIC")
champ.QC() # Alternatively QC.GUI(arraytype="EPIC")
myNorm <- champ.norm(beta=myLoad$beta, arraytype="EPIC",cores=30)
myNorm <- champ.norm(beta=myLoad$beta, method="PBC")

champ.SVD()
# If Batch detected, run champ.runCombat() here.This data is not suitable.
myDMP <- champ.DMP(arraytype="EPIC")

DMP.GUI()
myDMR <- champ.DMR()
DMR.GUI()
myDMR <- champ.DMR(arraytype="EPIC")
DMR.GUI(arraytype="EPIC")
myBlock <- champ.Block(arraytype="EPIC")
Block.GUI(arraytype="EPIC") # For this simulation data, not Differential Methylation Block is detected.
myGSEA <- champ.GSEA(arraytype="EPIC")

# champ.CNA(arraytype="EPIC")
# champ.CNA() function call for intensity data, which is not included in our Simulation data.





