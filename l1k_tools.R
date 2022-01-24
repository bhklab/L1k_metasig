library(fgsea)
library(utils)
library(R.utils)
library(coop)
library(cmapR)
library(data.table)
library(dplyr)


download_l1k_data <- function(mypath){
  # Download L1000 files from the GEO repository GSE92742
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FLevel5%5FCOMPZ%2EMODZ%5Fn473647x12328%2Egctx%2Egz", 
                destfile=paste(mypath, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz", sep="/"))
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FREADME%2Epdf", 
                destfile=paste(mypath, "GSE92742_Broad_LINCS_README.pdf", sep="/"))
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fcell%5Finfo%2Etxt%2Egz", 
                destfile=paste(mypath, "GSE92742_Broad_LINCS_cell_info.txt.gz", sep="/"))
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fgene%5Finfo%2Etxt%2Egz", 
                destfile=paste(mypath, "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="/"))
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fpert%5Finfo%2Etxt%2Egz", 
                destfile=paste(mypath, "GSE92742_Broad_LINCS_pert_info.txt.gz", sep="/"))
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fpert%5Fmetrics%2Etxt%2Egz", 
                destfile=paste(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt.gz", sep="/"))
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fsig%5Finfo%2Etxt%2Egz", 
                destfile=paste(mypath, "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="/"))
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fsig%5Fmetrics%2Etxt%2Egz", 
                destfile=paste(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt.gz", sep="/"))
  
  gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz", sep="/"))
  gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_cell_info.txt.gz", sep="/"))
  gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_gene_info.txt.gz", sep="/"))
  gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_pert_info.txt.gz", sep="/"))
  gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt.gz", sep="/"))
  gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_sig_info.txt.gz", sep="/"))
  gunzip(destfile=paste(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt.gz", sep="/"))
}


read_l1k_meta <- function(mypath, version=2017){
  if (version == 2017){
    mypath <- file.path(mypath, "2017")
    cellinfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_cell_info.txt", sep="/"), stringsAsFactors=FALSE)
    geneinfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_gene_info.txt", sep="/"), stringsAsFactors=FALSE)
    pertinfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_pert_info.txt", sep="/"), stringsAsFactors=FALSE)
    siginfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_sig_info.txt", sep="/"), stringsAsFactors=FALSE)
    instinfo <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_inst_info.txt", sep="/"), stringsAsFactors=FALSE)
    
    geneinfo$pr_gene_id <- as.character(geneinfo$pr_gene_id)
    
    pertmetrics <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_pert_metrics.txt", sep="/"), stringsAsFactors=FALSE)
    sigmetrics <- read.delim(file=paste(mypath, "GSE92742_Broad_LINCS_sig_metrics.txt", sep="/"), stringsAsFactors=FALSE)
    
    landmarks <- geneinfo[which(geneinfo$pr_is_lm == 1), ]
    
    return(list(cellinfo=cellinfo, 
                geneinfo=geneinfo,
                pertinfo=pertinfo,
                siginfo=siginfo,
                instinfo=instinfo,
                pertmetrics=pertmetrics,
                sigmetrics=sigmetrics, 
                landmarks=landmarks))
  } else if (version == 2020){
    mypath <- file.path(mypath, "2020")
    cellinfo <- read.delim(file=file.path(mypath, "cellinfo_beta.txt"), stringsAsFactors = FALSE)
    geneinfo <- read.delim(file=file.path(mypath, "geneinfo_beta.txt"), stringsAsFactors = FALSE)
    pertinfo <- read.delim(file=file.path(mypath, "compoundinfo_beta.txt"), stringsAsFactors = FALSE)
    siginfo <- read.delim(file=file.path(mypath, "siginfo_beta.txt"), stringsAsFactors = FALSE)
    
    geneinfo$gene_id <- as.character(geneinfo$gene_id)
    # Remap geneinfo - because the 2020 data is beta, the column names may change.  I map them to the same ones as 2017 for consistency within Leo:
    colnames(geneinfo)[match(c("gene_id", "gene_symbol", "gene_title"), colnames(geneinfo))] <- c("pr_gene_id", "pr_gene_symbol", "pr_gene_title")
    geneinfo$pr_is_lm <- as.numeric(geneinfo$feature_space == "landmark")
    geneinfo$pr_is_bing <- as.numeric(geneinfo$feature_space %in% c("landmark", "best inferred"))
    
    # Remap siginfo, cellinfo, pertinfo
    colnames(siginfo)[match(c("cmap_name", "cell_iname"), colnames(siginfo))] <- c("pert_iname", "cell_id")
    colnames(cellinfo)[match(c("cell_iname"), colnames(cellinfo))] <- c("cell_id")
    colnames(pertinfo)[match(c("cmap_name"), colnames(cellinfo))] <- c("pert_iname")
    
    landmarks <- geneinfo[which(geneinfo$feature_space == "landmark"), ]
    
    return(list(cellinfo=cellinfo, 
                geneinfo=geneinfo, 
                pertinfo=pertinfo, 
                siginfo=siginfo,
                landmarks=landmarks))
  }
}



# Compute distance is a function for computing similarities 
compute_distance <- function(ds, gene_annot, upgenes=c(), dngenes=c(), geneweights=0, 
                             geneset_name="", metric="cmap_score", 
                             gseaParam=1, nperms=1000){
  
  # Add condition to check input, upgenes, dngenes, or a weighted vector geneweights
  # Add condition to check metric, currently only cmap_score supported. 
  if (nperms <= 0){
    nperms <- 10
  }
  
  upgenes <- intersect(upgenes, gene_annot$pr_gene_symbol)
  dngenes <- intersect(dngenes, gene_annot$pr_gene_symbol)
  common <- intersect(upgenes, dngenes)
  
  upgenes <- setdiff(upgenes, common)
  dngenes <- setdiff(dngenes, common)
  
  upid <- geneinfo$pr_gene_id[match(upgenes, geneinfo$pr_gene_symbol)]
  dnid <- geneinfo$pr_gene_id[match(dngenes, geneinfo$pr_gene_symbol)]
  
  res <- data.frame(cid=character(length(ds@cid)), 
                    geneset=character(length(ds@cid)), 
                    pval=numeric(length(ds@cid)), 
                    padj=numeric(length(ds@cid)), 
                    ES=numeric(length(ds@cid)), 
                    NES=numeric(length(ds@cid)), 
                    size=numeric(length(ds@cid)), stringsAsFactors=FALSE)
  
  for (ii in seq(length(ds@cid))){
    if (ii %% 1000 == 0) {print(ii)}
    tempGSEA <- fgsea(list(upid=upid, dnid=dnid), ds@mat[, ii], nperm=nperms, gseaParam=gseaParam)
    combES <- (tempGSEA$ES[1] - tempGSEA$ES[2])/2 * abs(sign(tempGSEA$ES[1]) - sign(tempGSEA$ES[2]))/2
    combNES <- (tempGSEA$NES[1] - tempGSEA$NES[2])/2 * abs(sign(tempGSEA$NES[1]) - sign(tempGSEA$NES[2]))/2
    
    res[ii,] <- list(ds@cid[ii], geneset_name, 0, 0, combES, combNES, tempGSEA$size[1] + tempGSEA$size[2])
  }
  
  # Add a step to normalize with respect to touchstone set
  
  
  res
}

# This function should probably be deleted, I think it is legacy
compute_distance <- function(ds, gene_annot, siginfo, upgenes=c(), dngenes=c(), geneweights=0, 
                             geneset_name="mygeneset", metric="cmap_score", 
                             gseaParam=1, nperms=10, 
                             saveFiles=FALSE, savePath="."){
  
  # Add condition to check input, upgenes, dngenes, or a weighted vector geneweights
  # Add condition to check metric, currently only cmap_score supported. 
  if (nperms <= 0){
    nperms <- 1
  }
  
  upgenes <- intersect(upgenes, gene_annot$pr_gene_symbol)
  dngenes <- intersect(dngenes, gene_annot$pr_gene_symbol)
  common <- intersect(upgenes, dngenes)
  
  upgenes <- setdiff(upgenes, common)
  dngenes <- setdiff(dngenes, common)
  
  upid <- geneinfo$pr_gene_id[match(upgenes, geneinfo$pr_gene_symbol)]
  dnid <- geneinfo$pr_gene_id[match(dngenes, geneinfo$pr_gene_symbol)]
  
  if (metric == "cmap_score"){
    res <- data.frame(cid=character(length(ds@cid)), 
                      geneset=character(length(ds@cid)), 
                      pval=numeric(length(ds@cid)), 
                      padj=numeric(length(ds@cid)), 
                      ES=numeric(length(ds@cid)), 
                      NES=numeric(length(ds@cid)), 
                      GSEA_up=numeric(length(ds@cid)),
                      GSEA_dn=numeric(length(ds@cid)),
                      GSEA_up_p=numeric(length(ds@cid)),
                      GSEA_dn_p=numeric(length(ds@cid)),
                      size=numeric(length(ds@cid)), stringsAsFactors=FALSE)
    
    if (length(intersect(dnid, ds@rid)) > 0 & length(intersect(upid, ds@rid)) > 0){
      # Two-tailed query
      for (ii in seq(length(ds@cid))){
        if (ii %% 1000 == 0) {print(ii)}
        tempGSEA <- fgsea(list(upid=upid, dnid=dnid), ds@mat[, ii], nperm=nperms, gseaParam=gseaParam)
        
        combES <- (tempGSEA$ES[1] - tempGSEA$ES[2])/2 * abs(sign(tempGSEA$ES[1]) - sign(tempGSEA$ES[2]))/2
        combNES <- (tempGSEA$NES[1] - tempGSEA$NES[2])/2 * abs(sign(tempGSEA$NES[1]) - sign(tempGSEA$NES[2]))/2
        
        # This breaks for some reason for Altay with error:
        # Error: Values in 'p' vector must not contain NAs.
        # The correction below is supposed to correct this. 
        if (sum(is.na(tempGSEA$pval)) > 0){
          combp <- data.frame(p=1)
        } else {
          combp <- fisher(tempGSEA$pval)
        }
        
        res[ii,] <- list(ds@cid[ii], geneset_name, 
                         ifelse(combES == 0, 1, combp$p), 0, 
                         combES, combNES, 
                         tempGSEA$ES[1], tempGSEA$ES[2], 
                         tempGSEA$pval[1], tempGSEA$pval[2], tempGSEA$size[1] + tempGSEA$size[2])
      }
    } else if (length(intersect(upid, ds@rid)) > 0){
      # Single-tailed, up genes only
      for (ii in seq(length(ds@cid))){
        if (ii %% 1000 == 0) {print(ii)}
        tempGSEA <- fgsea(list(upid=upid), ds@mat[, ii], nperm=nperms, gseaParam=gseaParam)
        
        combES <- tempGSEA$ES[1]
        combNES <- tempGSEA$NES[1]
        combp <- tempGSEA$pval
        
        res[ii,] <- list(ds@cid[ii], geneset_name, 
                         ifelse(combES == 0, 1, combp), 0, 
                         combES, combNES, 
                         tempGSEA$ES[1], 0, 
                         tempGSEA$pval[1], tempGSEA$size[1])
      }
    } else if (length(intersect(dnid, ds@rid)) > 0){
      # Single-tailed, down genes only
      for (ii in seq(length(ds@cid))){
        if (ii %% 1000 == 0) {print(ii)}
        tempGSEA <- fgsea(list(dnid=dnid), ds@mat[, ii], nperm=nperms, gseaParam=gseaParam)
        
        combES <- -tempGSEA$ES[1]
        combNES <- -tempGSEA$NES[1]
        combp <- tempGSEA$pval
        
        res[ii,] <- list(ds@cid[ii], geneset_name, 
                         ifelse(combES == 0, 1, combp), 0, 
                         combES, combNES, 
                         tempGSEA$ES[1], 0, 
                         tempGSEA$pval[1], tempGSEA$size[1])
      }
    }
  } 
  else if (metric == "cosine"){
    res <- data.frame(cid=ds@cid, 
                      geneset=rep(geneset_name, length(ds@cid)),
                      metric=rep(metric, length(ds@cid)),
                      sim=numeric(length(ds@cid)))
    res$sim <- sapply(seq(dim(ds@mat)[2]), FUN=function(x) coop::cosine(ds@mat[,x], geneweights))
  }
  
  if (saveFiles){
    mydirname <- saveQuery(res, savePath, geneset_name, upgenes, dngenes, upid, dnid)
  }
  
  res
}


query_sig <- function(mysig_id, dspath, myrid, mycid, geneinfo, siginfo, savePath=".", tspath="", metric="cmap_score"){
  
  if (metric == "cmap_score"){
    mysets <- get_topk_genes(dspath, sig_id=mysig_id, geneinfo, k=50, geneset=myrid)
    ds <- parse_gctx(dspath, rid=myrid, cid = mycid)
    
    res <- compute_distance(ds, geneinfo, siginfo, 
                            upgenes=mysets$upset, 
                            dngenes=mysets$dnset,
                            geneweights=0, 
                            geneset_name=sprintf("%s_wtcs", mysig_id),
                            metric=metric,
                            saveFiles=TRUE, 
                            savePath=".")
    
  } else if (metric == "cosine"){
    
    ds_sig <- parse_gctx(dspath, cid=mysig_id, rid=myrid)
    gix <- which(colRanks(abs(ds_sig@mat)) > (dim(ds_sig@mat)[1] - 100))
    geneweights <- ds_sig@mat[gix]
    genenames <- geneinfo$pr_gene_symbol[match(ds_sig@rid[gix], geneinfo$pr_gene_id)]
    
    ds <- parse_gctx(dspath, rid=myrid, cid = mycid)
    
    
    res <- compute_distance(ds, geneinfo, siginfo, 
                            geneweights=geneweights, )
    
    res <- batch_query(dspath=dspath,
                       myrid=myrid,
                       mycid=mycid,
                       geneinfo, siginfo,
                       geneweights=geneweights,
                       genenames=genenames,
                       geneset_name=sprintf("%s_cos_n=100", mysig_id), 
                       metric=metric,
                       ts_sigs=ts_sigs, 
                       focus_sigs=focus_sigs,
                       saveFiles = TRUE, 
                       savePath = savePath,
                       tspath=tspath,
                       truncatequery=FALSE)
  }
}

compute_dist_matrix <- function(dspath, cid, rid, metric){
  
}


compute_cs_block <- function(ds1, ds2, k=50, gseaParam=1, nperms=1){
  # This function computes wtcs(x, y) where x is a signature drawn from ds1, and y one from ds2.
  # The up and down genesets from each element of ds1 are checked for enrichment in ds2. 
  upsets <- get_top_k(ds1@mat, k=k, decreasing=TRUE)
  dnsets <- get_top_k(ds1@mat, k=k, decreasing=FALSE)
  allsets <- c(upsets, dnsets)
  names(allsets) <- as.character(seq(length(allsets)))
  ncol <- dim(ds1@mat)[2]
  
  es <- matrix(numeric(dim(ds2@mat)[2] * ncol), ncol = ncol)
  printf("ds1 mat dimension: %d x %d \n", dim(ds1@mat)[1], dim(ds1@mat)[2])
  printf("ds2 mat dimension: %d x %d \n", dim(ds2@mat)[1], dim(ds2@mat)[2])
  
  for (ii in seq(dim(ds2@mat)[2])){
    tempGSEA <- fgsea(allsets, ds2@mat[, ii], nperm=nperms, gseaParam=gseaParam)
    # This generates p-values, but is more expensive. 
    # tGSEA <- fgseaMultilevel(allsets, ds2@mat[, ii], sampleSize=k)
    upscore <- tempGSEA$ES[1:ncol]
    dnscore <- tempGSEA$ES[(1+ncol):(2*ncol)]
    es[ii, ] <- (upscore - dnscore)/2 * abs(sign(upscore) - sign(dnscore))/2
    
    if (ii %% 100 == 0){
      printf("%d / %d\n", ii, dim(ds2@mat)[2])
    }
  }
  
  printf("ES Names dimension: %d x %d \n", length(ds1@cid), length(ds2@cid))
  rownames(es) <- ds2@cid
  colnames(es) <- ds1@cid
  
  printf("Returning ES, class = %s \n", class(es))
  return(es)
}

compute_cosine_block <- function(ds1, ds2, k=0){
  # This function computes cosine(x, y), where x is a signature drawn from ds1 and y fromone from ds2. 
  # It returns an n1 x n2 matrix, where n1 is the number of columns in ds1, and n2 is the number in ds2.
  # Cosine distance is symmetric.
  # This can probably be optimized - or at least parallelized. 
  
  if (k > 0){
    ds1@mat <- as.matrix(ds1@mat * (colRanks(ds1@mat, preserveShape = TRUE) <= k | colRanks(ds1@mat, preserveShape = TRUE) > (dim(ds1@mat)[1] - k)))
  }
  
  cosres <- matrix(numeric(dim(ds1@mat)[2] * dim(ds2@mat)[2]), nrow = dim(ds1@mat)[2], dimnames=list(ds1@cid, ds2@cid))
  # Currently iterative; parallelize
  for (ii in seq(dim(ds1@mat)[2])){
    cosres[ii,] <- sapply(seq(dim(ds2@mat)[2]), FUN=function(x) cosine(ds1@mat[,ii], ds2@mat[,x]))
  }
  return(cosres)
}


get_level5_ds <- function(mypath){
  f <- list.files(mypath, pattern="level5.*gctx", ignore.case=TRUE)
  
  if (length(f) != 1){
    stop(sprintf("In get_level5_ds: Number of matched files in %s is %d, but should be 1", mypath, length(f)))
  } 
  return(file.path(mypath, f[1]))
}


get_pert_sigs <- function(dspath, siginfo, sigmetrics, geneinfo, pert_id="", pert_iname="", 
                          cell_ids=list(), saveFiles=FALSE, outpath="."){
  # Identify signatures of interest
  if (pert_iname != ""){
    mysigs <- siginfo[grep(pert_iname, siginfo$pert_iname, ignore.case=TRUE),]
  } else {
    mysigs <- siginfo[grep(pert_id, siginfo$pert_id, ignore.case=TRUE),]
  }
  
  if (length(cell_ids) != 0){
    # Restrict cell lines of interest
    mysigs <- mysigs[which(mysigs$cell_id %in% cell_ids), ]
  }
  
  mysigm <- sigmetrics[match(mysigs$sig_id, sigmetrics$sig_id),]
  mysigs <- cbind(mysigs, mysigm[, setdiff(colnames(mysigm), colnames(mysigs))])
  
  if (dim(mysigs)[1] == 0){
    print(sprintf("Error: No signatures found for pert_id = %s, pert_iname = %s", pert_id, pert_iname))
    return(list())
  }
  
  ds <- parse_gctx(dspath, cid=mysigs$sig_id)
  dsmat <- ds@mat
  rownames(dsmat) <- geneinfo$pr_gene_symbol[match(rownames(dsmat), geneinfo$pr_gene_id)]
  
  # Do some averaging or something
  ds_avg <- rowMeans(dsmat)
  ds_best <- dsmat[, which(mysigs$distil_cc_q75 == max(mysigs$distil_cc_q75))[1]]
  ds_strong <- NA
  if (sum(mysigs$distil_cc_q75 > 0.2) > 0){
    ds_strong <- rowMeans(dsmat[, which(mysigs$distil_cc_q75 > 0.2), drop=FALSE])
  }
  ds_bycell <- sapply(unique(mysigs$cell_id), FUN = function(x) rowMeans(dsmat[ ,which(mysigs$cell_id == x), drop=FALSE]))
  
  agg_sigs <- cbind("average"=ds_avg, "best"=ds_best, "reproducible avg"=ds_strong, ds_bycell)
  
  # This function seems quite slow, despite its simplicity.  
  #ds_bycell <- t(as.data.table(t(dsmat)) %>% group_by(mysigs$cell_id) %>% summarise_all(mean))
  
  # Dose response ??
  mysigs <- mysigs[, !(colnames(mysigs) %in% c("distil_id"))]
  
  if (saveFiles){
    write.table(x = dsmat, paste(outpath, sprintf("sigs_%s%s.txt", pert_id, pert_iname), sep="/"), sep="\t")
    write.table(x = mysigs, paste(outpath, sprintf("siginfo_%s%s.txt", pert_id, pert_iname), sep="/"), sep="\t")
    write.table(x = agg_sigs, paste(outpath, sprintf("sigs_aggregate_%s%s.txt", pert_id, pert_iname), sep="/"), sep="\t")
  }
  
  return(list(sigs=dsmat, aggregate_sigs=agg_sigs, siginfo=mysigs))
}


get_gene_perts <- function(dspath, siginfo, geneinfo, sigmetrics, gene_id="", cell_ids=list(), saveFiles=FALSE, outpath=".", pert_type=c("trt_cp")){
  geneix <- match(gene_id, geneinfo$pr_gene_symbol)
  if (is.na(geneix)){
    print(sprintf("Error: no gene found for gene_id = %s.  Check spelling.", gene_id))
    return(0)
  } 
  
  print(geneinfo[geneix,])
  # To do: Print quality of gene if inferred.
  
  sigs <- siginfo$sig_id[siginfo$pert_type %in% pert_type]
  ds <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneix], cid=sigs)
  
  expr_mat <- t(as.data.table(ds@mat))
  
  sigstr <- sigmetrics[match(rownames(expr_mat), sigmetrics$sig_id), c("sig_id", "distil_ss")]
  expr_mat <- cbind(expr_mat, expr_mat/sigstr$distil_ss)
  
  colnames(expr_mat) <- c(gene_id, sprintf("%s_norm", gene_id))
  expr_mat <- cbind(expr_mat, siginfo[match(rownames(expr_mat), siginfo$sig_id), c("sig_id", "pert_id", "pert_iname", "cell_id", "pert_idose", "pert_itime")], 
                    sigmetrics[match(rownames(expr_mat), sigmetrics$sig_id), c("distil_cc_q75", "distil_ss")])
  
  summary_mat <- as.data.frame(expr_mat %>% group_by(pert_iname) %>% dplyr::summarise(mean_expr=mean(get(gene_id)), mean_norm=mean(get(sprintf("%s_norm", gene_id)))))
  
  expr_mat <- expr_mat[order(expr_mat[,sprintf("%s_norm", gene_id)], decreasing=TRUE),]
  summary_mat <- summary_mat[order(summary_mat[, "mean_norm"], decreasing=TRUE),]
  
  if (saveFiles){
    write.table(expr_mat, file=paste(outpath, sprintf("%s_exprvalues_n=%d.txt", gene_id, dim(expr_mat)[1]), sep="/"), sep="/t")
    write.table(summary_mat, file=paste(outpath, sprintf("%s_exprbypert_n=%d.txt", gene_id, dim(summary_mat)[1]), sep="/"), sep="/t")
  }
  
  return(list(expr_mat=expr_mat, summary_mat=summary_mat))
}

compute_cs_mats <- function(mat1, mat2, kgenes=50, gseaParam=1, nperms=1, parallel=0, numCores=detectCores()){
  if (parallel == 0){
    upsets <- get_top_k(mat1, k=kgenes, decreasing=TRUE)
    dnsets <- get_top_k(mat1, k=kgenes, decreasing=FALSE)
    allsets <- c(upsets, dnsets)
    names(allsets) <- as.character(seq(length(allsets)))
    ncol <- dim(mat1)[2]
    
    es <- matrix(numeric(dim(mat2)[2] * ncol), ncol = ncol)
    
    for (ii in seq(dim(mat2)[2])){
      tempGSEA <- fgsea(allsets, mat2[, ii], nperm=nperms, gseaParam=gseaParam)
      # This generates p-values, but is more expensive. 
      # tGSEA <- fgseaMultilevel(allsets, ds2@mat[, ii], sampleSize=k)
      upscore <- tempGSEA$ES[1:ncol]
      dnscore <- tempGSEA$ES[(1+ncol):(2*ncol)]
      es[ii, ] <- (upscore - dnscore)/2 * abs(sign(upscore) - sign(dnscore))/2
    }
    
    rownames(es) <- colnames(mat2)
    colnames(es) <- colnames(mat1)
    return(es)
    
  } else {
    upsets <- get_top_k(mat1, k=kgenes, decreasing=TRUE)
    dnsets <- get_top_k(mat1, k=kgenes, decreasing=FALSE)
    allsets <- c(upsets, dnsets)
    names(allsets) <- as.character(seq(length(allsets)))
    ncol <- dim(mat1)[2]
    
    tempGSEA <- mclapply(seq(dim(mat2)[2]), FUN=function(x) fgsea(allsets, mat2[,x], nperm=nperms, gseaParam=gseaParam), mc.cores=numCores)
    zes <- lapply(tempGSEA, FUN=function(y) (y$ES[1:ncol] - y$ES[(1+ncol):(2*ncol)])/2 * abs(sign(y$ES[1:ncol]) - sign(y$ES[(1+ncol):(2*ncol)]))/2)
    es <- t(matrix(unlist(zes), ncol = dim(mat2)[2]))
    
    rownames(es) <- colnames(mat2)
    colnames(es) <- colnames(mat1)
    return(es)
  }
}


get_top_k <- function(dsmat, k=50, decreasing=TRUE){
  
  if (k > 500){
    k <- 500
  }
  mysets <- lapply(seq(dim(dsmat)[2]), FUN=function(x) rownames(dsmat)[order(dsmat[,x], decreasing=decreasing)[seq(k)]])
  return(mysets)
}


compute_sim_block <- function(ds1, ds2, metric="cosine", kgenes=0, gseaParam=1, nperms=1, parallel=0, numCores=detectCores()){
  
  if (!all.equal(ds1@rid, ds2@rid)){
    return("Error: input datasets did not have the same row space or rids.")
  }
  
  if (parallel == 0){
    
    if (metric == "cosine"){
      
      if (kgenes > 0){
        ds1@mat <- as.matrix(ds1@mat * (colRanks(ds1@mat, preserveShape = TRUE) <= kgenes | colRanks(ds1@mat, preserveShape = TRUE) > (dim(ds1@mat)[1] - kgenes)))
      }
      
      cosres <- matrix(numeric(dim(ds1@mat)[2] * dim(ds2@mat)[2]), nrow = dim(ds1@mat)[2], dimnames=list(ds1@cid, ds2@cid))
      # Currently iterative; parallelize
      for (ii in seq(dim(ds1@mat)[2])){
        cosres[ii,] <- sapply(seq(dim(ds2@mat)[2]), FUN=function(x) cosine(ds1@mat[,ii], ds2@mat[,x]))
      }
      cosres <- t(cosres)
      colnames(cosres) <- ds1@cid
      rownames(cosres) <- ds2@cid
      return(cosres)
      
    } else if (metric == "wtcs"){
      
      upsets <- get_top_k(ds1@mat, k=kgenes, decreasing=TRUE)
      dnsets <- get_top_k(ds1@mat, k=kgenes, decreasing=FALSE)
      allsets <- c(upsets, dnsets)
      names(allsets) <- as.character(seq(length(allsets)))
      ncol <- dim(ds1@mat)[2]
      
      es <- matrix(numeric(dim(ds2@mat)[2] * ncol), ncol = ncol)
      #printf("ds1 mat dimension: %d x %d \n", dim(ds1@mat)[1], dim(ds1@mat)[2])
      #printf("ds2 mat dimension: %d x %d \n", dim(ds2@mat)[1], dim(ds2@mat)[2])
      
      for (ii in seq(dim(ds2@mat)[2])){
        tempGSEA <- fgsea(allsets, ds2@mat[, ii], nperm=nperms, gseaParam=gseaParam)
        # This generates p-values, but is more expensive. 
        # tGSEA <- fgseaMultilevel(allsets, ds2@mat[, ii], sampleSize=k)
        upscore <- tempGSEA$ES[1:ncol]
        dnscore <- tempGSEA$ES[(1+ncol):(2*ncol)]
        es[ii, ] <- (upscore - dnscore)/2 * abs(sign(upscore) - sign(dnscore))/2
        
        if (ii %% 100 == 0){
          printf("%d / %d\n", ii, dim(ds2@mat)[2])
        }
      }
      
      rownames(es) <- ds2@cid
      colnames(es) <- ds1@cid
      return(es)
      
    } else if (metric == "pearson") {
      
      pearson_res <- cor(ds1@mat, ds2@mat, method="pearson")
      pearson_res <- t(pearson_res)
      colnames(pearson_res) <- ds1@cid
      rownames(pearson_res) <- ds2@cid
      return(pearson_res)
      
    } else if (metric == "spearman") {
      
      spearman_res <- cor(ds1@mat, ds2@mat, method="spearman")
      spearman_res <- t(spearman_res)
      colnames(spearman_res) <- ds1@cid
      rownames(spearman_res) <- ds2@cid
      return(spearman_res)
      
    }
  } else {
    # Parallelized
    
    if (metric == "cosine"){
      
      if (kgenes > 0){
        ds1@mat <- as.matrix(ds1@mat * (colRanks(ds1@mat, preserveShape = TRUE) <= kgenes | colRanks(ds1@mat, preserveShape = TRUE) > (dim(ds1@mat)[1] - kgenes)))
      }
      
      # Currently iterative; parallelize
      
      x <- mclapply(seq(dim(ds1@mat)[2]), FUN=function(x) {
        sapply(seq(dim(ds2@mat)[2]), FUN=function(y) cosine(ds1@mat[,x], ds2@mat[,y]))
      }, mc.cores=numCores)
      cosres <- matrix(unlist(x), nrow=dim(ds1@mat)[2])
      
      cosres <- t(cosres)
      colnames(cosres) <- ds1@cid
      rownames(cosres) <- ds2@cid
      
      return(cosres)
    }  else if (metric == "wtcs"){
      
      upsets <- get_top_k(ds1@mat, k=kgenes, decreasing=TRUE)
      dnsets <- get_top_k(ds1@mat, k=kgenes, decreasing=FALSE)
      allsets <- c(upsets, dnsets)
      names(allsets) <- as.character(seq(length(allsets)))
      ncol <- dim(ds1@mat)[2]
      
      #printf("ds1 mat dimension: %d x %d \n", dim(ds1@mat)[1], dim(ds1@mat)[2])
      #printf("ds2 mat dimension: %d x %d \n", dim(ds2@mat)[1], dim(ds2@mat)[2])
      
      tempGSEA <- mclapply(seq(dim(ds2@mat)[2]), FUN=function(x) fgsea(allsets, ds2@mat[,x], nperm=nperms, gseaParam=gseaParam), mc.cores=numCores)
      zes <- lapply(tempGSEA, FUN=function(y) (y$ES[1:ncol] - y$ES[(1+ncol):(2*ncol)])/2 * abs(sign(y$ES[1:ncol]) - sign(y$ES[(1+ncol):(2*ncol)]))/2)
      es <- t(matrix(unlist(zes), ncol = dim(ds2@mat)[2]))
      
      rownames(es) <- ds2@cid
      colnames(es) <- ds1@cid
      return(es)
      
    } else if (metric == "pearson") {
      
      pearson_res <- cor(ds2@mat, ds1@mat, method="pearson")
      colnames(pearson_res) <- ds1@cid
      rownames(pearson_res) <- ds2@cid
      return(pearson_res)
      
    } else if (metric == "spearman") {
      
      spearman_res <- cor(ds2@mat, ds1@mat, method="spearman")
      colnames(spearman_res) <- ds1@cid
      rownames(spearman_res) <- ds2@cid
      return(spearman_res)
      
    }
  }
}