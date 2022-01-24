library(cmapR)
library(ggplot2)
library(parallel)
library(png)

datapath <- "~/Work/bhk/data/l1k"
outdir <- "~/Work/bhk/analysis/l1k/202107_highd/"
dspath <- file.path(datapath, "2017", "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")

#### Load Metadata
l1kmeta <- read_l1k_meta(datapath, version=2017)
attach(l1kmeta)


#### Frankensig correlation vs sig size
pertcount5 <- table(siginfo$pert_iname[siginfo$pert_type == "trt_cp"])
myperts <- names(sort(pertcount5, decreasing=TRUE)[1:200])

pertdf <- data.frame()
frankendf <- list()

dsnull <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=sample(siginfo$sig_id[siginfo$pert_type == "trt_cp"], 20000))
dsdmso <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=sample(siginfo$sig_id[siginfo$pert_type == "ctl_vehicle"], 10000))
franknull <- get_frankencorr(dsnull@mat, nmax=400, iter=1000)
frankdmso <- get_frankencorr(dsdmso@mat, nmax=400, iter=1000)

for (mypert in myperts){
  print(sprintf("%s: %d/%d", mypert, match(mypert, myperts), length(myperts))) 
  
  mysigs <- siginfo$sig_id[siginfo$pert_iname == mypert]
  myccq75 <- sigmetrics[sigmetrics$sig_id %in% mysigs, "distil_cc_q75"]
  pertccq75 <- mean(myccq75[myccq75 > -1])
  
  ds <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=mysigs)
  frankpert <- get_frankencorr(ds@mat, nmax=400, iter=100)
  summarydf <- rbind(cbind(frankpert$cordf, data=mypert), cbind(franknull$cordf, data="Random Cpds"), cbind(frankdmso$cordf, data="DMSO"))
  
  #pdf(file=file.path(outdir, "frankenplots", sprintf("%s_%d.pdf", mypert, length(mysigs))), width=10, height=8)
  #print(ggplot(summarydf, aes(x=frankensize, y=meancor, group=data, color=data)) + geom_line() + 
  #  geom_errorbar(aes(ymin=meancor-sdcor, ymax=meancor+sdcor), width=0.2) + theme_minimal() + ylim(-0.2, 1) + 
  #  labs(title=sprintf("Correlation of Split Average Signatures for %s, N=%d; mean ccq75=%0.2f", mypert, length(mysigs), pertccq75), 
  #       x="Number of signatures in average", y="Mean Pearson Correlation"))
  #dev.off()
  
  frankendf <- c(frankendf, list(frankpert$cordf))
  
  myl <- max(length(frankpert$cordf$meancor), length(rownames(pertdf)))
  pertdf <- rbind(pertdf, head(c(frankpert$cordf$meancor, NA*seq_along(myl)), myl))
  if (is.null(colnames(pertdf))){
    colnames(pertdf) <- frankpert$cordf$frankensize
  }
}

rownames(pertdf) <- myperts
colnames(pertdf) <- frankendf[[1]]$frankensize
names(frankendf) <- myperts

#saveRDS(file=file.path(outdir, "frankensummary", "pertdf_pearson.RDS"), pertdf)
#saveRDS(file=file.path(outdir, "frankensummary", "frankendf_pearson.RDS"), frankendf)

saveRDS(file=file.path(outdir, "frankensummary", "pertdf_pearson_n=400.rds"), pertdf)
saveRDS(file=file.path(outdir, "frankensummary", "frankendf_pearson_n=400.rds"), frankendf)


# For the paper figure:

dsnull <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=sample(siginfo$sig_id[siginfo$pert_type == "trt_cp"], 10000))
dsdmso <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=sample(siginfo$sig_id[siginfo$pert_type == "ctl_vehicle"], 10000))
franknull <- get_frankencorr(dsnull@mat, nmax=100, iter=1000)
frankdmso <- get_frankencorr(dsdmso@mat, nmax=100, iter=1000)


#### Frankensig correlation vs sig size - other metrics
pertcount5 <- table(siginfo$pert_iname[siginfo$pert_type == "trt_cp"])
myperts <- names(sort(pertcount5, decreasing=TRUE)[1:200])

for (mymetric in c("wtcs")){ #c("spearman", "cosine", "wtcs")){
  simdf <- data.frame(pert_iname=character(), n=numeric(), ccq75=numeric(), sim5=numeric(), sim10=numeric(), sim25=numeric())
  
  iter <- 100 #ifelse(mymetric == "wtcs", 25, 100)
  for (mypert in myperts){
    print(sprintf("%s %s: %d/%d", mymetric, mypert, match(mypert, myperts), length(myperts))) 
    
    mysigs <- siginfo$sig_id[siginfo$pert_iname == mypert]
    myccq75 <- sigmetrics[sigmetrics$sig_id %in% mysigs, "distil_cc_q75"]
    pertccq75 <- mean(myccq75[myccq75 > -1])
    
    ds <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=mysigs)
    dsnull <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=sample(siginfo$sig_id[siginfo$pert_type == "trt_cp"], 1000))
    dsdmso <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=sample(siginfo$sig_id[siginfo$pert_type == "ctl_vehicle"], 1000))
    
    frankpert <- get_frankensim(ds@mat, nmax=55, iter=iter, metric=mymetric)
    franknull <- get_frankensim(dsnull@mat, nmax=55, iter=iter, metric=mymetric)
    frankdmso <- get_frankensim(dsdmso@mat, nmax=55, iter=iter, metric=mymetric)
    
    summarydf <- rbind(cbind(frankpert$cordf, data=mypert), cbind(franknull$cordf, data="Random Cpds"), cbind(frankdmso$cordf, data="DMSO"))
    resdf <- rbind(cbind(reshape2::melt(frankpert$simmat), dataset=mypert), 
                   cbind(reshape2::melt(franknull$simmat), dataset="Random Cpds"), 
                   cbind(reshape2::melt(frankdmso$simmat), dataset="DMSO"))
    
    if (metric == "wtcs"){
      pdf(file=file.path(outdir, sprintf("frankenplots_%s", mymetric), sprintf("%s_fullbox_%d.pdf", mypert, length(mysigs))), width=10, height=8)
      print(ggplot(resdf, aes(x=as.factor(Var1), y=value, fill=dataset, color=dataset)) + geom_jitter(position=position_jitterdodge()) + 
              geom_boxplot(alpha=0.2, compact=TRUE) + geom_dotplot() + theme_minimal() + ylim(-0.2, 1) + 
              labs(title=sprintf("%s Similarity of Split Average Signatures for %s, N=%d; mean ccq75=%0.2f", mymetric, mypert, length(mysigs), pertccq75), 
                   x="Number of signatures in average", y=sprintf("Mean %s", mymetric)))
      dev.off()
    } 
    
    pdf(file=file.path(outdir, sprintf("frankenplots_%s", mymetric), sprintf("%s_%d.pdf", mypert, length(mysigs))), width=10, height=8)
    print(ggplot(summarydf, aes(x=frankensize, y=meancor, group=data, color=data)) + geom_line() + 
            geom_errorbar(aes(ymin=meancor-sdcor, ymax=meancor+sdcor), width=0.2) + theme_minimal() + ylim(-0.2, 1) + 
            labs(title=sprintf("%s Similarity of Split Average Signatures for %s, N=%d; mean ccq75=%0.2f", mymetric, mypert, length(mysigs), pertccq75), 
                 x="Number of signatures in average", y=sprintf("Mean %s", mymetric)))
    dev.off()
    
    simdf <- rbind(simdf, data.frame(pert_iname=mypert, n=length(mysigs), ccq75=pertccq75, 
                                       sim5=frankpert$cordf$meancor[1], sim10=frankpert$cordf$meancor[2], sim25=frankpert$cordf$meancor[5])) 
  }
  
  saveRDS(file=file.path(outdir, "frankensummary", sprintf("pertdf_%s.RDS", mymetric)), simdf)
}




# Load pertdfs:
frankendf <- readRDS(file=file.path(outdir, "frankensummary", "frankendf_pearson.RDS"))
pertdf <- readRDS(file=file.path(outdir, "frankensummary", "pertdf.RDS"))



#### Megapert df

topperts <- names(sort(pertcount5, decreasing=TRUE)[1:20])

for (mypert in topperts){
  print(sprintf("%s: %d/%d", mypert, match(mypert, topperts), length(topperts))) 
  
  mysigs <- siginfo$sig_id[siginfo$pert_iname == mypert]
  myccq75 <- sigmetrics[sigmetrics$sig_id %in% mysigs, "distil_cc_q75"]
  pertccq75 <- mean(myccq75[myccq75 > -1])
  
  ds <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=mysigs)
  dsnull <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=sample(siginfo$sig_id[siginfo$pert_type == "trt_cp"], 2000))
  dsdmso <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=sample(siginfo$sig_id[siginfo$pert_type == "ctl_vehicle"], 2000))
  
  frankpert <- get_frankencorr(ds@mat, nmax=500, iter=100)
  franknull <- get_frankencorr(dsnull@mat, nmax=500, iter=100)
  frankdmso <- get_frankencorr(dsdmso@mat, nmax=500, iter=100)
  
  summarydf <- rbind(cbind(frankpert$cordf, data=mypert), cbind(franknull$cordf, data="Random Cpds"), cbind(frankdmso$cordf, data="DMSO"))
  
  pdf(file=file.path(outdir, "frankenplots", sprintf("megaplot_%s_%d.pdf", mypert, length(mysigs))), width=10, height=8)
  print(ggplot(summarydf, aes(x=frankensize, y=meancor, group=data, color=data)) + geom_line() + 
          geom_errorbar(aes(ymin=meancor-sdcor, ymax=meancor+sdcor), width=0.2) + theme_minimal() + ylim(-0.2, 1) + 
          labs(title=sprintf("Correlation of Split Average Signatures for %s, N=%d; mean ccq75=%0.2f", mypert, length(mysigs), pertccq75), 
               x="Number of signatures in average", y="Mean Pearson Correlation"))
  dev.off()
}


#### Do the k-frankensigs self-correlate?

pertcount5 <- table(siginfo$pert_iname[siginfo$pert_type == "trt_cp"])
myperts <- names(pertcount5[pertcount5 >= 60])

nvals <- c(1, 3, 5, 10, 20, 30)
subperts <- sample(myperts, 30)
mycormats <- list()

for (n in nvals){
  print(n)
  pertdata <- list()
  
  for (mypert in subperts){
    print(mypert)
    mysigs <- siginfo$sig_id[siginfo$pert_iname == mypert]
    ds <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=mysigs)
    frankpert <- get_frankensigs(ds@mat, k=n, returnk=1)
    
    pertdata <- c(pertdata, list(frankpert))
  }
  
  names(pertdata) <- subperts
  
  cormat <- matrix(numeric(length(pertdata) * length(pertdata)), nrow=length(pertdata))
  sdmat <- matrix(numeric(length(pertdata) * length(pertdata)), nrow=length(pertdata))
  
  for (ii in seq_along(pertdata)){
    mycorrs <- cor(pertdata[[ii]])[upper.tri(cor(pertdata[[ii]]))]
    cormat[ii,ii] <- mean(mycorrs)
    sdmat[ii,ii] <- sd(mycorrs)
    
    for (jj in ii + seq_len(length(pertdata) - ii)){
      mycorrs <- cor(pertdata[[ii]], pertdata[[jj]])
      cormat[ii,jj] <- mean(mycorrs)
      sdmat[ii,jj] <- sd(mycorrs)
    }
  }
  
  mycormats <- c(mycormats, list(list(cormat=cormat, sdmat=sdmat)))
}

pdf(file.path(outdir, "frankensummary", "summary_crossfranken_boxplot_n=5-30.pdf"), width=10, height=8)
boxplot(diag(mycormats[[1]]$cormat), mycormats[[1]]$cormat[upper.tri(mycormats[[1]]$cormat)], 
        diag(mycormats[[2]]$cormat), mycormats[[2]]$cormat[upper.tri(mycormats[[2]]$cormat)], 
        diag(mycormats[[3]]$cormat), mycormats[[3]]$cormat[upper.tri(mycormats[[3]]$cormat)], 
        diag(mycormats[[4]]$cormat), mycormats[[4]]$cormat[upper.tri(mycormats[[4]]$cormat)], 
        names=c("Self n=5", "Cross n=5", "Self n=10", "Cross n=10", "Self n=20", "Cross n=20", "Self n=30", "Cross n=30"),
        col=c("red", "blue"), ylab="Pearson Correlation", main="Correlation of Frankensignatures within and across perturbation")
grid()
dev.off()

pdf(file.path(outdir, "frankensummary", "summary_franken10_vs_ccq75.pdf"), width=10, height=8)
plot(pertdf$ccq75, pertdf$pearson10, pch=16, xlab="Average Replicate correlation", ylab="Mean Franken10 Correlation", 
     main="Frankencorrelation vs average replicate correlation")
dev.off()



####  Do maximal frankensigs correlate?

myperts <- names(sort(pertcount5, decreasing=TRUE)[1:1000])

pertcounts <- sort(pertcount5, decreasing=TRUE)[1:1000]
fperts <- matrix(numeric(length(myperts)*978), nrow=978, dimnames=list(geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], myperts))
fcorrs <- list(length(myperts))

# FIX: Change frankpert from full frankensigs to half frankensigs to match fcorrs
for (ii in seq_along(myperts)){ 
  mypert <- myperts[ii]
  print(mypert)
  mysigs <- siginfo$sig_id[siginfo$pert_iname == mypert]
  ds <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=mysigs)
  
  frankpert <- get_frankensigs(ds@mat, k=(dim(ds@mat)[2]), returnk=0, return2=0)
  fcorrs[ii] <- list(bootstrap_maxfrankencorr(ds@mat, iter=100))
  fperts[,ii] <- frankpert
}

mycors <- cor(fperts)

pdf(file.path(outdir, "frankensummary", "fullfrankensigs_corr_densityplot.pdf"), width=10, height=8)
plot(density(sapply(fcorrs, mean), bw=0.01), col="red", xlim=c(-1,1), xlab="Pearson Correlation", 
     main="Self correlation vs cross correlation for Top 1000 compound Frankensignatures", lwd=2)
lines(density(mycors[upper.tri(mycors)], bw=0.01), col="blue", lwd=2)
legend(x="topleft", legend=c("Same pert, Franken 2-split correlation", "Cross franken correlation, different cpd"), col=c("red", "blue"), lwd=c(2,2))
grid()
dev.off()



#### Compute general similarity matrices  ####
csigs <- siginfo$sig_id[siginfo$pert_type == "trt_cp"]

ds <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=csigs)

mymetrics <- c("cosine", "pearson", "spearman", "wtcs")
for (met in mymetrics){
  print(met)
  print(system.time(compute_sim_block(subset_gct(ds, cid=1:100), subset_gct(ds, cid=1:100), metric=met, parallel=1)))
}

# csigs_sample <- csigs[sample(length(csigs), 500)]
simds <- readRDS(file.path(outdir, "simmats", "cosine_lm_n=205034x500.rds"))
csigs_sample <- colnames(simds)

for (met in mymetrics){
  print(met)
  mysim <- compute_sim_block(subset_gct(ds, cid=csigs_sample), ds, metric=met, parallel=1)
  saveRDS(mysim, file = file.path(outdir, "simmats", sprintf("%s_lm_n=%dx%d.rds", met, length(csigs), length(csigs_sample))))
}

# Memory efficient computation of wtcs similarity block, if the above doesn't work:
mysim <- numeric()
tic()
for (ii in seq(1, 206)){
  print(ii)
  tsim <- compute_sim_block(subset_gct(ds, cid=csigs_sample), subset_gct(ds, cid=ds@cid[((ii-1)*1000+1):(min(length(ds@cid), ii*1000))]), 
                            metric="wtcs", parallel=1, kgenes=50)
  mysim <- rbind(mysim, tsim)
}
toc()
saveRDS(mysim, file=file.path(outdir, "simmats", sprintf("wtcs_lm_n=%dx%d.rds", dim(mysim)[1], dim(mysim)[2])))



#### Null metasignature correlations ####
csigs <- siginfo$sig_id[siginfo$pert_type == "trt_cp"]
ds <- parse_gctx(dspath, rid=geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], cid=csigs)

pertcount5 <- table(siginfo$pert_iname[siginfo$pert_type == "trt_cp"])

sig_agg <- sigmetrics[sigmetrics$pert_type == "trt_cp", ] %>% group_by(pert_iname) %>% 
  summarize(nsample = sum(distil_nsample), mycount = length(distil_nsample), pert_type = max(pert_type))

pertcorr <- data.frame(pert_iname=names(pertcount5), sig_count=as.numeric(pertcount5), 
                       inst_count=sig_agg$nsample[match(names(pertcount5), sig_agg$pert_iname)], corr=numeric(length(myperts)))
myperts <- pertcorr$pert_iname

avgsig <- rowMeans(ds@mat)
for (ii in seq_along(myperts)){
  if (ii %% 100 == 0){print(ii)}
  mypert <- myperts[ii]
  pertsigs <- siginfo$sig_id[siginfo$pert_iname == mypert]
  if (length(pertsigs) > 1){
    mysig <- rowMeans(ds@mat[, match(pertsigs, ds@cid)])
  } else {
    mysig <- ds@mat[, match(pertsigs, ds@cid)]
  }
  othersig <- (avgsig * length(ds@cid) - mysig * pertcorr$sig_count[ii]) / (length(ds@cid) - pertcorr$sig_count[ii])
  pertcorr$corr[ii] <- cor(mysig, othersig)
}

saveRDS(pertcorr, file=file.path(outdir, "cor_anomaly", "metasig_pert_corr.rds"))

############################
##### Manuscript plots #####
############################

pertcorr <- readRDS(file=file.path(outdir, "cor_anomaly", "metasig_pert_corr.rds"))

# Figure Example Autocorrelation: Tozasertib, median compound
# Get the franken20 correlation:
f20 <- sapply(frankendf, FUN=function(x) x$meancor[x$frankensize == 20])

# Choose the 10th, 101st, and 190th rankeds compounds to be representative. 
# These were chosen to roughly correspond to 5th, 50th, and 95th quantiles, and drugs
# with alphanumeric pert_inames were excluded in favor of recognition. 
# Compounds: Tolazamide (blood glucose lowering drug), Tozasertib (AURK inhibitor), and Vorinostat (HDAC inhibitor)

figperts <- c("estradiol", "tozasertib", "vorinostat")
figdf <- rbind(cbind(franknull$cordf, name="All Cpds"), 
               cbind(frankdmso$cordf, name="DMSO"), 
               cbind(frankendf[[figperts[1]]], name=figperts[1]),
               cbind(frankendf[[figperts[2]]], name=figperts[2]),
               cbind(frankendf[[figperts[3]]], name=figperts[3]))

pdf(file.path(outdir, "figures/excpd_frankencorr_pearson.pdf"), width=8, height=7)
ggplot(figdf, aes(x=frankensize, y=meancor, group=name)) + geom_line(aes(colour=factor(name)), size=1.5) + theme_minimal() +
  geom_ribbon(aes(ymin=meancor - sdcor, ymax=meancor + sdcor, fill=factor(name)), alpha=0.3) +
  scale_fill_manual(guide="none", name="", values=c("vorinostat"="#E69F00", "tozasertib"="#56B4E9", "estradiol"="#009E73", "All Cpds"="black", "DMSO"="brown")) +
  scale_colour_manual(name="", values=c("vorinostat"="#E69F00", "tozasertib"="#56B4E9", "estradiol"="#009E73", "All Cpds"="black", "DMSO"="brown")) + 
  xlab("K averaged signatures") + ylab("PCC") + 
  theme(legend.text=element_text(size=14), text=element_text(size=16), legend.position=c(0.87,0.25), 
        legend.background=element_rect(fill="white", size=0, linetype="blank"))
dev.off()

# Figure: all compound autocorrelation 
pdf(file.path(outdir, "figures/fig1b_allcpd_frankencorr_pearson.pdf"), width=8, height=6)
matplot(c(5, 10, 25), t(pertdf[, 4:6]), type="l", col="grey", xlab="K, number of averaged signatures", ylab="Pearson Correlation", main="Increasing autocorrelation of K-signatures for 200 compounds", ylim=c(0,1))
points(c(5,10,25), franknull$cordf$meancor[c(5, 10, 25)], col="forestgreen", lty=2, type="b")
points(c(5,10,25), frankdmso$cordf$meancor[c(5, 10, 25)], col="red", lty=2, type="b")
dev.off()

pdf(file.path(outdir, "figures/allcpd_frankencorr_pearson_full.pdf"), width=8, height=7)
par(mar=c(3.5,3.2,1,1))
plot(frankendf[[1]]$frankensize, frankendf[[1]]$meancor, col="#0072B2", type="l", xlab="", ylab="", ylim=c(0,1), lwd=1.5)
#, xlab="K, number of averaged signatures", ylab="Pearson Correlation", main="", ylim=c(0,1))
title(ylab="PCC", xlab="K, number of averaged signatures", line=2, cex.lab=1.4)
for(ii in seq(2, length(frankendf))){
  lines(frankendf[[ii]]$frankensize, frankendf[[ii]]$meancor, col="#0072B2", type="l", lwd=1.5)
}
lines(franknull$cordf$frankensize, franknull$cordf$meancor, col="#000000", type="l", lwd=2.5, lty=6)
lines(frankdmso$cordf$frankensize, frankdmso$cordf$meancor, col="#D55E00", type="l", lwd=2.5, lty=2)
legend(x="bottomright", legend=c("Compound-specific", "Random compounds", "DMSOs"), 
       col=c("#0072B2", "#000000", "#D55E00"), lty=c(1,6,2), lwd=c(2,2,2))
dev.off()


# Big all compound autocorrelation
pdf(file.path(outdir, "figures/allcpd_frankencorr_pearson_full_x=400_log.pdf"), width=8, height=7)
plot(frankendf[[1]]$frankensize, frankendf[[1]]$meancor, col="#0072B2", type="l", xlab="", ylab="", ylim=c(0,1), lwd=1.5, log="x")
#, xlab="K, number of averaged signatures", ylab="Pearson Correlation", main="", ylim=c(0,1))
title(ylab="PCC", xlab="K, number of averaged signatures", line=2, cex.lab=1.4)
for(ii in seq(2, length(frankendf))){
  lines(frankendf[[ii]]$frankensize, frankendf[[ii]]$meancor, col="#0072B2", type="l", lwd=1.5)
}
lines(franknull$cordf$frankensize, franknull$cordf$meancor, col="#000000", type="l", lwd=2.5, lty=6)
lines(frankdmso$cordf$frankensize, frankdmso$cordf$meancor, col="#D55E00", type="l", lwd=2.5, lty=2)
legend(x="bottomright", legend=c("Compound-specific", "Random compounds", "DMSOs"), 
       col=c("#0072B2", "#000000", "#D55E00"), lty=c(1,6,2), lwd=c(2,2,2))
dev.off()

# Figure: metric similarity
coslm <- readRDS(file.path(outdir, "simmats/cosine_lm_n=205034x500.rds"))
pearlm <- readRDS(file.path(outdir, "simmats/pearson_lm_n=205034x500.rds"))
sprmlm <- readRDS(file.path(outdir, "simmats/spearman_lm_n=205034x500.rds"))
wtcslm <- readRDS(file.path(outdir, "simmats/wtcs_lm_n=205034x500.rds"))

rk_coslm <- apply(coslm, 2, rank)/dim(coslm)[1]
rk_pearlm <- apply(pearlm, 2, rank)/dim(pearlm)[1]
rk_sprmlm <- apply(sprmlm, 2, rank)/dim(sprmlm)[1]
rk_wtcslm <- apply(wtcslm, 2, rank)/dim(wtcslm)[1]

mysigm <- sigmetrics[match(rownames(coslm), sigmetrics$sig_id), ]
mynsample <- sapply(mysigm$distil_nsample, FUN=function(x) ceiling(min(x, 10)/2))
grnsm <- sapply(mysigm$distil_nsample, FUN=function(x) ceiling(min(x, 10)))
binsm <- numeric(length(grnsm))
binsm[grnsm %in% c(1,2)] <- 1
binsm[grnsm == 3] <- 2
binsm[grnsm >= 4] <- 3

mycols <- palette.colors(5)
breaks <- seq(-1, 1, 0.01)
pdf(file.path(outdir, "figures/metric_nsample_cosine_cdf.pdf"), width=10, height=8)
plot(breaks, ecdf_pointwise(coslm[mynsample == 1,], breaks=breaks)[,2], col=mycols[1], xlab="Cosine", 
     main="Cosine Distance", type="l", lwd=2, xlim=c(-0.5, 0.5), ylab="Cumulative Density")
for (ii in seq(2, 5)) {
  lines(breaks, ecdf_pointwise(coslm[mynsample == ii,], breaks=breaks)[,2], col=mycols[ii], lwd=2)
}
legend(x="bottomright", legend=c("Nsample = 1-2", "3-4", "5-6", "7-8", "9+"), lwd=2, col=mycols)
dev.off()

pdf(file.path(outdir, "figures/metric_nsample_pearson_cdf.pdf"), width=10, height=8)
plot(breaks, ecdf_pointwise(pearlm[mynsample == 1,], breaks=breaks)[,2], col=mycols[1], xlab="Pearson", 
     main="Pearson", type="l", lwd=2, xlim=c(-0.5, 0.5), ylab="Cumulative Density")
for (ii in seq(2, 5)) {
  lines(breaks, ecdf_pointwise(pearlm[mynsample == ii,], breaks=breaks)[,2], col=mycols[ii], lwd=2)
}
legend(x="bottomright", legend=c("Nsample = 1-2", "3-4", "5-6", "7-8", "9+"), lwd=2, col=mycols)
dev.off()

pdf(file.path(outdir, "figures/metric_nsample_spearman_cdf.pdf"), width=10, height=8)
plot(breaks, ecdf_pointwise(sprmlm[mynsample == 1,], breaks=breaks)[,2], col=mycols[1], xlab="Spearman", 
     main="Spearman", type="l", lwd=2, xlim=c(-0.5, 0.5), ylab="Cumulative Density")
for (ii in seq(2, 5)) {
  lines(breaks, ecdf_pointwise(sprmlm[mynsample == ii,], breaks=breaks)[,2], col=mycols[ii], lwd=2)
}
legend(x="bottomright", legend=c("Nsample = 1-2", "3-4", "5-6", "7-8", "9+"), lwd=2, col=mycols)
dev.off()

pdf(file.path(outdir, "figures/metric_nsample_wtcs_cdf.pdf"), width=10, height=8)
plot(breaks, ecdf_pointwise(wtcslm[which(mynsample == 1),], breaks=breaks)[,2], col=mycols[1], xlab="WTCS", 
     main="WTCS", type="l", lwd=2, xlim=c(-1, 1), ylab="Cumulative Density")
for (ii in seq(2, 5)) {
  lines(breaks, ecdf_pointwise(wtcslm[which(mynsample == ii),], breaks=breaks)[,2], col=mycols[ii], lwd=2)
}
legend(x="bottomright", legend=c("Nsample = 1-2", "3-4", "5-6", "7-8", "9+"), lwd=2, col=mycols)
dev.off()

# Use binned nsamples:
pdf(file.path(outdir, "figures/metric_nsample_cosine_bin_cdf.pdf"), width=10, height=8)
plot(breaks, ecdf_pointwise(coslm[binsm == 1,], breaks=breaks)[,2], col=mycols[1], xlab="Cosine", 
     main="Cosine Distance", type="l", lwd=2, xlim=c(-0.4, 0.4), ylab="Cumulative Density")
for (ii in seq(2, 3)) {
  lines(breaks, ecdf_pointwise(coslm[binsm == ii,], breaks=breaks)[,2], col=mycols[ii], lwd=2)
}
legend(x="bottomright", legend=c("Nsample = 1-2", "3", "4+"), lwd=2, col=mycols)
dev.off()

pdf(file.path(outdir, "figures/metric_nsample_pearson_bin_cdf.pdf"), width=10, height=8)
plot(breaks, ecdf_pointwise(pearlm[binsm == 1,], breaks=breaks)[,2], col=mycols[1], xlab="Pearson", 
     main="Pearson", type="l", lwd=2, xlim=c(0.1, 0.4), ylim=c(0.8, 1), ylab="Cumulative Density")
for (ii in seq(2, 3)) {
  lines(breaks, ecdf_pointwise(pearlm[binsm == ii,], breaks=breaks)[,2], col=mycols[ii], lwd=2)
}
legend(x="bottomright", legend=c("Nsample = 1-2", "3", "4+"), lwd=2, col=mycols)
dev.off()

pdf(file.path(outdir, "figures/metric_nsample_spearman_bin_cdf.pdf"), width=10, height=8)
plot(breaks, ecdf_pointwise(sprmlm[binsm == 1,], breaks=breaks)[,2], col=mycols[1], xlab="Spearman", 
     main="Spearman", type="l", lwd=2, xlim=c(-0.4, 0.4), ylab="Cumulative Density")
for (ii in seq(2, 3)) {
  lines(breaks, ecdf_pointwise(sprmlm[binsm == ii,], breaks=breaks)[,2], col=mycols[ii], lwd=2)
}
legend(x="bottomright", legend=c("Nsample = 1-2", "3", "4+"), lwd=2, col=mycols)
dev.off()

pdf(file.path(outdir, "figures/metric_nsample_wtcs_bin_cdf_inset.pdf"), width=8, height=6)
par(mar=c(3,3,1,1))
plot(breaks, ecdf_pointwise(wtcslm[binsm[1:(dim(wtcslm)[1])] == 1,], breaks=breaks)[,2], col=mycols[1], xlab="", 
     main="", type="l", lwd=2.5, xlim=c(0.2, 0.6), ylim=c(0.5, 1), ylab="")
title(ylab="Cumulative Density", xlab="WTCS", line=1.9, cex.lab=1.4)
for (ii in seq(2, 3)) {
  lines(breaks, ecdf_pointwise(wtcslm[binsm[1:(dim(wtcslm)[1])] == ii,], breaks=breaks)[,2], col=mycols[ii], lwd=2.5)
}
legend(x="bottomright", legend=c("Nsample = 1-2", "3", "4+"), lwd=2.5, col=mycols)
dev.off()


# Figure: proportion in the top 1%
query1pct <- rbind(
  data.frame(metric="cosine", nsample=c("1-2", "3", "4+"), rate=sapply(seq(3), 
      FUN=function(x) mean(rk_coslm[binsm == x,] >= 0.99))/1e-2),
  data.frame(metric="pearson", nsample=c("1-2", "3", "4+"), rate=sapply(seq(3), 
      FUN=function(x) mean(rk_pearlm[binsm == x,] >= 0.99))/1e-2),
  data.frame(metric="spearman", nsample=c("1-2", "3", "4+"), rate=sapply(seq(3), 
      FUN=function(x) mean(rk_sprmlm[binsm == x,] >= 0.99))/1e-2),
  data.frame(metric="wtcs", nsample=c("1-2", "3", "4+"), rate=sapply(seq(3), 
      FUN=function(x) mean(rk_wtcslm[binsm == x,] >= 0.99))/1e-2)
)

pdf(file.path(outdir, "figures/metric_querybias_bar_1pct.pdf"), width=8, height=6)
ggplot(query1pct, aes(fill=nsample, y=rate, x=metric)) + geom_bar(position="dodge", stat="identity") + 
  theme_minimal() + xlab("") + ylab("Frequency") + 
  theme(text=element_text(size=18))
dev.off()

wtcs_test12 <- wilcox.test(wtcslm[binsm == 1, ], wtcslm[binsm == 2, ])
wtcs_test13 <- wilcox.test(wtcslm[binsm == 1, ], wtcslm[binsm == 3, ])
wtcs_test23 <- wilcox.test(wtcslm[binsm == 2, ], wtcslm[binsm == 3, ])


# Figure: Distribution of number of signatures per compound
pdf(file.path(outdir, "figures/fig1d_compound_nsample_hist.pdf"), width=8, height=6)
hist(sapply(pertcount5, FUN=function(x) min(x, 50)), breaks=50, xlab="Nsample", ylab="Frequency", main="")
dev.off()

# Figure 1E: 
# pdf(file.path(outdir, "figures/fig1e_replicatecorr_vs_frankencorr.pdf"), width=8, height=6)
# plot(pertdf$ccq75, pertdf$pearson10, pch=16, xlab="ccq75 - replicate correlation", ylab="Mean 10-metasignature Pearson", main="")
# dev.off()


# Meta signature correlation:
# Overlay DLEPS Figure 2i:
myimg <- readPNG(file.path(outdir, "figures/dleps_fig2i.png"))
pertcorr$highsigs <- pertcorr$sig_count >= 15

jx <- which(pertcorr$inst_count > 5)
ix <- which(pertcorr$sig_count >= 15)
kx <- which(pertcorr$inst_count > 5 & pertcorr$sig_count < 15)
hidens <- density(pertcorr$corr[ix], bw=0.01)
lodens <- density(pertcorr$corr[kx], bw=0.01)

pdf(file.path(outdir, "figures/metasig_cor.pdf"), width=8, height=6)
par(mar=c(3,3,1,1))
plot(c(-0.55, 1), c(0, 2.5), type="n", main="", xlab="", ylab="")
title(ylab="Density", xlab="PCC", line=1.9, cex.lab=1.4)
rasterImage(myimg, xleft = -0.58, xright=1.052, ybottom=0, ytop=2.5)
lines(density(pertcorr$corr[jx], bw=0.01), lwd=3, col="darkgrey")
lines(hidens$x, length(ix)/length(jx)*hidens$y, lwd=3, col="magenta")
lines(lodens$x, length(kx)/length(jx)*lodens$y, lwd=3, col="black")
legend(x="topleft", legend=c("All compounds", "High cpds", "Low cpds", "DLEPS Fig 2i"), col=c("grey", "magenta", "black", "blue"), lwd=c(3,3,3,5))
dev.off()

pdf(file.path(outdir, "figures/metasig_cor_vs_sig.pdf"), width=8, height=6)
par(mar=c(3,3,1,1))
ggplot(pertcorr[jx,], aes(x=log10(sig_count), y=corr)) + geom_point(aes(color=highsigs)) + geom_smooth() + 
  scale_color_manual(values=c("black", "magenta")) + theme_minimal() + 
  xlab("Signature Count") + ylab("PCC") + scale_y_continuous(breaks=seq(-0.3, 0.9, 0.3), labels=seq(-0.3, 0.9, 0.3), limits=c(-0.3, 0.9)) +
  theme(text=element_text(size=18), legend.position="none") + 
  scale_x_continuous(breaks=c(log10(3), 1, log10(30), 2, log10(300), 3), 
                     labels=c(3, 10, 30, 100, 300, 1000))
dev.off()



