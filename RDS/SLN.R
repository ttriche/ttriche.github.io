library(SingleCellExperiment) # this is needed to load the data
library(compositions) # this is needed for centering ADT logratios

skip_reload <- TRUE
if (skip_reload == TRUE) { 

  # load the data off of GitHub (pre-smallified, split, and merged)
  SLN <- readRDS(url("https://ttriche.github.io/RDS/SLN_merged.rds"))

} else { # if we want to go through all the steps 

  # load the data off of GitHub (pre-smallified but not merged)

  # assays(SLN_111)$X <- as(assays(SLN_111)$X, "CsparseMatrix")
  SLN_111 <- readRDS(url("https://ttriche.github.io/RDS/SLN_111_SCE.rds"))
  assayNames(SLN_111) <- "counts"
  sapply(assays(SLN_111), class)
  #      counts 
  # "dgCMatrix" 

  # assays(SLN_206)$X <- as(assays(SLN_206)$X, "CsparseMatrix")
  SLN_206 <- readRDS(url("https://ttriche.github.io/RDS/SLN_206_SCE.rds"))
  assayNames(SLN_206) <- "counts"
  sapply(assays(SLN_206), class)
  #      counts 
  # "dgCMatrix" 


  # label the antibodies
  antibodies <- list(
    SLN_111=metadata(SLN_111)$protein_names,
    SLN_206=metadata(SLN_206)$protein_names
  )

  # rename the cells if any of them clash (a few dozen cells do)
  if (length(intersect(colnames(SLN_111), colnames(SLN_206))) > 0) { 
    colnames(SLN_111) <- paste0("SLN_111_", colnames(SLN_111))
    rownames(reducedDims(SLN_111)$protein_expression) <- colnames(SLN_111)
    colnames(SLN_206) <- paste0("SLN_206_", colnames(SLN_206))
    rownames(reducedDims(SLN_111)$protein_expression) <- colnames(SLN_111)
  }

  # label the counts for each of these antibodies in each dataset
  colnames(reducedDims(SLN_111)$protein_expression) <- antibodies$SLN_111
  colnames(reducedDims(SLN_206)$protein_expression) <- antibodies$SLN_206

  # note the HTO (hash-tag oligo) counts: 
  lapply(antibodies, function(x) grep("HTO", setdiff(merged_ab, x), val=TRUE))
  # 
  # $SLN_111
  # [1] "HTO_B6_LN_r4_206_A0301"  "HTO_B6_spl_r4_206_A0302"
  #
  # $SLN_206
  # [1] "HTO_B6_spl_r4_206_A0301" "HTO_B6_LN_r4_206_A0302" 
  
  # merge the reducedDim slots holding antibody-derived tag counts, move em out
  SLN_111_tags <- split.data.frame(t(reducedDim(SLN_111, "protein_expression")),
                                   substr(antibodies$SLN_111, 1, 3)) # ADT / HTO
  SLN_111_tags <- lapply(SLN_111_tags, as, "CsparseMatrix") # like the mains 
  # bonus: way easier to sanity-check!
  # Note how SLN_111_tags$HTO corresponds to SLN_111$hash_id
  #
  # Lymph Node   Negative     Spleen 
  #       7606         15       9207 

  # merge the reducedDim slots holding antibody-derived tag counts, move em out
  SLN_206_tags <- split.data.frame(t(reducedDim(SLN_206, "protein_expression")),
                                   substr(antibodies$SLN_206, 1, 3)) # ADT / HTO
  SLN_206_tags <- lapply(SLN_206_tags, as, "CsparseMatrix") # like the mains 
  # bonus: way easier to sanity-check!
  # Note how SLN_206_tags$HTO corresponds to SLN_206$hash_id
  table(SLN_206$hash_id)
  #
  #  Doublet Lymph Node   Negative     Spleen 
  #       24       7368         29       8399 

  ADT_names <- c(rownames(SLN_111_tags$ADT), rownames(SLN_206_tags$ADT))
  merged_ADTs <- unique(ADT_names)

  # SLN_206 ADTs are a strict superset of those in SLN_111
  length(setdiff(merged_ADTs, rownames(SLN_206_tags$ADT)))
  # [1] 0
  length(setdiff(merged_ADTs, rownames(SLN_111_tags$ADT)))
  # [1] 97

  # label the sources 
  SLN_111$source <- "SLN_111"
  SLN_206$source <- "SLN_206"

  # drop the reducedDims and merge the mRNA UMI counts
  # (reducedDims is not where the ADTs should go anyhow!)
  reducedDims(SLN_111) <- NULL 
  reducedDims(SLN_206) <- NULL 
  SLN <- cbind(SLN_111, SLN_206) 

  # drop doublets and negatives
  SLN <- SLN[, SLN$hash_id %in% c("Lymph Node", "Spleen")] 
  ADTs <- list(SLN_111=SLN_111_tags$ADT, SLN_206=SLN_206_tags$ADT)
  keep_good_ADTs <- function(x, SLN) x[, which(colnames(x) %in% colnames(SLN))]
  ADTs <- lapply(ADTs, keep_good_ADTs, SLN=SLN) 

  sapply(ADTs, ncol)
  # SLN_111 SLN_206 
  #   16813   15767 

  table(SLN$source)
  # SLN_111 SLN_206 
  #   16813   15767

  # turn the ADT counts into centered log-ratios (compositions::clr)
  library(compositions)
  
  # Seurat's CLR is not a true CLR, fwiw: 
  clr_function <- function(x) {
    log1p(x=x/(exp(x=sum(log1p(x=x[x > 0]), na.rm=TRUE)/length(x=x))))
  }

  # a proper CLR is log(x) - mean(log(x)), but this can at least be adapted
  cl1pr <- function(x) {
    log1p(x=x/(expm1(x=sum(log1p(x=x[x > 0]), na.rm=TRUE)/length(x=x))))
  }

  # vectorized for ADTs with rows as antibodies and columns as cells: 
  sCLR <- function(ADTs, how=c("cl1pr", "clr_function", "clr")) {
    as(apply(ADTs, 2, match.fun(how)), "CsparseMatrix")
  }

  # inverse
  sCLRinv <- function(z, how=c("cl1pr","clr_function","clr"), orig=gsi.orig(z)){
    how <- match.arg(how) 
    if (how == "clr") {
      return(apply(z, 2, clrInv))
    } else { 
      sInv <- function(w) acomp(gsi.recodeC2M(expm1(w), ninf=0, nan=NaN, na=NA))
      return(apply(z, 2, sInv))
    } 
  }

  # let's plot the results to see how this affects things 
  library(ComplexHeatmap) 
  CLRs <- c(correctedPos="cl1pr", Seurat="clr_function", compositional="clr")
  SLN_111_CLRs <- lapply(CLRs, sCLR, ADTs=ADTs$SLN_111)

  # let's look at dendritic cells 
  DCs <- colnames(SLN)[which(grepl("DC", SLN$cell_types))] 
  HMs <- list() 
  for (i in names(CLRs)) {
    cells <- intersect(colnames(SLN_111_CLRs[[i]]), DCs)
    clrmat <- as.matrix(SLN_111_CLRs[[i]][, cells])
    rownames(clrmat) <- sapply(strsplit(rownames(clrmat), "_"), `[`, 2)
    HMs[[i]] <- Heatmap(clrmat, name=paste(i, "CLR"),
                        show_column_names=FALSE, 
                        column_title=paste(i, "CLR"), 
                        row_names_gp=gpar(fontsize=7))
  }
  HMs[[1]] + HMs[[2]] + HMs[[3]]
  dev.copy2pdf(file="centeredLogRatioTransforms.pdf")
  # given that it is invertible, the left one seems like the sensible one

  # add NA rows for SLN_111 ADTs, merge, and prep for cl1pr 
  ADT <- as.matrix(ADTs$SLN_206)
  SLN_111_ADT_raw <- as.matrix(ADTs$SLN_111) 
  missing_111 <- setdiff(rownames(ADT), rownames(SLN_111_ADT_raw))
  SLN_111_NAs <- matrix(NA,
                        ncol=ncol(SLN_111_ADT_raw),
                        nrow=length(missing_111))
  rownames(SLN_111_NAs) <- missing_111
  colnames(SLN_111_NAs) <- colnames(SLN_111_ADT_raw)
  SLN_111_with_NAs <- rbind(SLN_111_ADT_raw, SLN_111_NAs)[rownames(ADT),] 
  merged_ADT <- cbind(SLN_111_with_NAs, ADT)[, colnames(SLN)]

  # add ADTs as an altExp
  altExp(SLN, "ADT") <- 
    SummarizedExperiment(assays=list(counts=as(merged_ADT, "CsparseMatrix"),
                                     CLR=as(sCLR(merged_ADT, how="cl1pr"), 
                                            "CsparseMatrix")),
                         colData=colData(SLN))
  mainExpName(SLN) <- "RNA"

  # label ADTs with their antigens
  rowData(altExp(SLN))$antigen <- 
    sapply(strsplit(rownames(altExp(SLN)), "(\\(|_)"), `[`, 2)

  # try and match up with gene names  
  library(Mus.musculus)
  rowData(altExp(SLN))$gene <- 
    mapIds(Mus.musculus, rowData(altExp(SLN))$antigen, "SYMBOL", "ALIAS")
  rowData(altExp(SLN))$hasRNA <- 
    !is.na(match(rowData(altExp(SLN))$gene, rownames(SLN)))

  # map RNA back to ADTs 
  rowData(SLN)$ADT <-
    rownames(altExp(SLN))[match(rownames(SLN), rowData(altExp(SLN))$gene)]
  rowData(SLN)$hasADT <- 
    !is.na(match(rownames(SLN), rowData(altExp(SLN))$gene)) 

  # the beefy one that GitHub won't take 
  saveRDS(SLN, file="SLN_merged_full.rds") 

  # downsample a bit so that GitHub *will* take it 
  # use harmony to merge the batches and 
  # label by cluster and cell type to preserve both
  library(velocessor)
  library(harmony)
  SLN_full <- logNormCounts(SLN)
  SLN_full$sample <- SLN_full$source
  SLN_full <- harmonize_velo_txis(SLN_full, "sample")
  colLabels(SLN_full) <- SLN_full$leiden_subclusters
  SLN_full$cellsourceclust <- paste0(SLN_full$source, ",",
                                     SLN_full$cell_types, ",",
                                     SLN_full$leiden_subclusters)
  saveRDS(SLN_full, file="SLN_merged_full.rds") 
 
  min(sort(table(SLN_full$cellsourceclust)))
  # 19
  
  samplecells <- function(obj, clusters, ideal=500) {
  
    pops <- sort(table(clusters))
    clusts <- names(pops) 
    graball <- names(which(pops <= ideal)) 
    keepall <- which(clusts %in% graball)

    grabsome <- setdiff(clusts, graball)
    tosample <- which(clusters %in% grabsome)
    samplesets <- split(tosample, clusters[tosample])
  
    keep <- c(keepall, do.call(c, lapply(samplesets, sample, size=ideal)))
    pct <- round((length(keep) / length(clusters))*100, 1)
    message("Kept ", length(keep), " (", pct, "%) of ",
            length(clusters), " cells in ", length(clusts), " clusters.")
    obj[, keep]

  }  
  SLN <- samplecells(SLN_full, SLN_full$cellsourceclust, 500) 

  # save the result and upload to github
  saveRDS(SLN, file="SLN_merged.rds")

}
