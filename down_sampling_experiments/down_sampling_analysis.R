library("Seurat")
library("EnImpute")
# source("EnImpute.R")
# source("alra.R")
source("downsamp_evaluation.R")

############################################################a#################
##  Impute the baron data
baron = readRDS("baron.rds")
baron_imputation_result = EnImpute(baron$count.samp)
saveRDS(baron_imputation_result, "baron_imputation_result.rds")
#baron_imputation_result = readRDS("baron_imputation_result.rds")

# Calculate correlation with reference data, CMD between the Pearson correlation 
# matrix derived from the imputed (and observed) matrix and the Pearson correlation 
# matrix derived from the reference matrix
K = length(baron_imputation_result$Methods.used)
ncell = dim(baron$count.samp)[2]
ngene = dim(baron$count.samp)[1]
baron_cor.gene = matrix(0,ngene,K+2)
baron_cor.cell = matrix(0,ncell,K+2)
baron_CMD.gene = rep(0,K+2)
baron_CMD.cell = rep(0,K+2)

ref.count =  log(normalizeData(baron$count.ref)+1)
ref.cor.gene = cor(t(ref.count))
ref.cor.cell = cor(ref.count)
for (k in 1:K){
  cat("k=",k)
  cat("\n")
  imputed.count = log(normalizeData(baron_imputation_result$count.imputed.individual[,,k])+1)
  baron_cor.gene[,k] = t(get.cor.gene(ref.count, imputed.count))
  baron_cor.cell[,k] = t(get.cor.cell(ref.count, imputed.count))
  imputed.cor.gene = cor(t(imputed.count))
  imputed.cor.gene[is.na(imputed.cor.gene)] = 0
  imputed.cor.cell = cor(imputed.count)
  imputed.cor.cell[is.na(imputed.cor.cell)] = 0
  baron_CMD.gene[k] = calc_cmd(ref.cor.gene, imputed.cor.gene)
  baron_CMD.cell[k] = calc_cmd(ref.cor.cell, imputed.cor.cell)
}

EnImpute.count =  log(normalizeData(baron_imputation_result$count.EnImpute.exp)+1)
baron_cor.gene[,K+1] = t(get.cor.gene(ref.count, EnImpute.count))
baron_cor.cell[,K+1] = t(get.cor.cell(ref.count, EnImpute.count))
EnImpute.cor.gene = cor(t(EnImpute.count))
EnImpute.cor.cell = cor(EnImpute.count)
baron_CMD.gene[K+1] = calc_cmd(ref.cor.gene, EnImpute.cor.gene)
baron_CMD.cell[K+1] = calc_cmd(ref.cor.cell, EnImpute.cor.cell)


sample.count =  log(normalizeData(baron$count.samp)+1)
baron_cor.gene[,K+2] = t(get.cor.gene(ref.count, sample.count))
baron_cor.cell[,K+2] = t(get.cor.cell(ref.count, sample.count))
sample.cor.gene = cor(t(sample.count))
sample.cor.cell = cor(sample.count)
baron_CMD.gene[K+2] = calc_cmd(ref.cor.gene, sample.cor.gene)
baron_CMD.cell[K+2] = calc_cmd(ref.cor.cell, sample.cor.cell)

save(baron, baron_imputation_result, baron_cor.gene, baron_cor.cell, baron_CMD.gene, baron_CMD.cell, file = "baron_result.RData")





############################################################a#################
## Impute the baron manno data
manno = readRDS("manno.rds")
manno_imputation_result = EnImpute(manno$count.samp)
saveRDS(manno_imputation_result, "manno_imputation_result.rds")
#manno_imputation_result = readRDS("manno_imputation_result.rds")

# Calculate correlation with reference data, CMD between the Pearson correlation 
# matrix derived from the imputed (and observed) matrix and the Pearson correlation 
# matrix derived from the reference matrix

K = length(manno_imputation_result$Methods.used)
ncell = dim(manno$count.samp)[2]
ngene = dim(manno$count.samp)[1]
manno_cor.gene = matrix(0,ngene,K+2)
manno_cor.cell = matrix(0,ncell,K+2)
manno_CMD.gene = rep(0,K+2)
manno_CMD.cell = rep(0,K+2)

ref.count =  log(normalizeData(manno$count.ref)+1)
ref.cor.gene = cor(t(ref.count))
ref.cor.cell = cor(ref.count)
for (k in 1:K){
  cat("k=",k)
  cat("\n")
  imputed.count = log(normalizeData(manno_imputation_result$count.imputed.individual[,,k])+1)
  manno_cor.gene[,k] = t(get.cor.gene(ref.count, imputed.count))
  manno_cor.cell[,k] = t(get.cor.cell(ref.count, imputed.count))
  imputed.cor.gene = cor(t(imputed.count))
  imputed.cor.gene[is.na(imputed.cor.gene)] = 0
  imputed.cor.cell = cor(imputed.count)
  imputed.cor.cell[is.na(imputed.cor.cell)] = 0
  manno_CMD.gene[k] = calc_cmd(ref.cor.gene, imputed.cor.gene)
  manno_CMD.cell[k] = calc_cmd(ref.cor.cell, imputed.cor.cell)
}

EnImpute.count =  log(normalizeData(manno_imputation_result$count.EnImpute.exp)+1)
manno_cor.gene[,K+1] = t(get.cor.gene(ref.count, EnImpute.count))
manno_cor.cell[,K+1] = t(get.cor.cell(ref.count, EnImpute.count))
EnImpute.cor.gene = cor(t(EnImpute.count))
EnImpute.cor.cell = cor(EnImpute.count)
manno_CMD.gene[K+1] = calc_cmd(ref.cor.gene, EnImpute.cor.gene)
manno_CMD.cell[K+1] = calc_cmd(ref.cor.cell, EnImpute.cor.cell)


sample.count =  log(normalizeData(manno$count.samp)+1)
manno_cor.gene[,K+2] = t(get.cor.gene(ref.count, sample.count))
manno_cor.cell[,K+2] = t(get.cor.cell(ref.count, sample.count))
sample.cor.gene = cor(t(sample.count))
sample.cor.cell = cor(sample.count)
manno_CMD.gene[K+2] = calc_cmd(ref.cor.gene, sample.cor.gene)
manno_CMD.cell[K+2] = calc_cmd(ref.cor.cell, sample.cor.cell)

save(manno, manno_imputation_result, manno_cor.gene, manno_cor.cell, manno_CMD.gene, manno_CMD.cell, file = "manno_result.RData")


############################################################a#################
## Impute the chen data
chen = readRDS("chen.rds")
chen_imputation_result = EnImpute(chen$count.samp, DrImpute = F)
saveRDS(chen_imputation_result, "chen_imputation_result.rds")
#chen_imputation_result = readRDS("chen_imputation_result.rds")

# Calculate correlation with reference data, CMD between the Pearson correlation 
# matrix derived from the imputed (and observed) matrix and the Pearson correlation 
# matrix derived from the reference matrix
K = length(chen_imputation_result$Methods.used)
ncell = dim(chen$count.samp)[2]
ngene = dim(chen$count.samp)[1]
chen_cor.gene = matrix(0,ngene,K+2)
chen_cor.cell = matrix(0,ncell,K+2)
chen_CMD.gene = rep(0,K+2)
chen_CMD.cell = rep(0,K+2)

ref.count =  log(normalizeData(chen$count.ref)+1)
ref.cor.gene = cor(t(ref.count))
ref.cor.cell = cor(ref.count)
for (k in 1:K){
  cat("k=",k)
  cat("\n")
  imputed.count = log(normalizeData(chen_imputation_result$count.imputed.individual[,,k])+1)
  chen_cor.gene[,k] = t(get.cor.gene(ref.count, imputed.count))
  chen_cor.cell[,k] = t(get.cor.cell(ref.count, imputed.count))
  imputed.cor.gene = cor(t(imputed.count))
  imputed.cor.gene[is.na(imputed.cor.gene)] = 0
  imputed.cor.cell = cor(imputed.count)
  imputed.cor.cell[is.na(imputed.cor.cell)] = 0
  chen_CMD.gene[k] = calc_cmd(ref.cor.gene, imputed.cor.gene)
  chen_CMD.cell[k] = calc_cmd(ref.cor.cell, imputed.cor.cell)
}

EnImpute.count =  log(normalizeData(chen_imputation_result$count.EnImpute.exp)+1)
chen_cor.gene[,K+1] = t(get.cor.gene(ref.count, EnImpute.count))
chen_cor.cell[,K+1] = t(get.cor.cell(ref.count, EnImpute.count))
EnImpute.cor.gene = cor(t(EnImpute.count))
EnImpute.cor.cell = cor(EnImpute.count)
chen_CMD.gene[K+1] = calc_cmd(ref.cor.gene, EnImpute.cor.gene)
chen_CMD.cell[K+1] = calc_cmd(ref.cor.cell, EnImpute.cor.cell)


sample.count =  log(normalizeData(chen$count.samp)+1)
chen_cor.gene[,K+2] = t(get.cor.gene(ref.count, sample.count))
chen_cor.cell[,K+2] = t(get.cor.cell(ref.count, sample.count))
sample.cor.gene = cor(t(sample.count))
sample.cor.cell = cor(sample.count)
chen_CMD.gene[K+2] = calc_cmd(ref.cor.gene, sample.cor.gene)
chen_CMD.cell[K+2] = calc_cmd(ref.cor.cell, sample.cor.cell)

save(chen, chen_imputation_result, chen_cor.gene, chen_cor.cell, chen_CMD.gene, chen_CMD.cell, file = "chen_result.RData")



############################################################a#################
## ##  Impute the zeisel data 
zeisel = readRDS("zeisel.rds")
zeisel_imputation_result = EnImpute(zeisel$count.samp)
saveRDS(zeisel_imputation_result, "zeisel_imputation_result.rds")
#zeisel_imputation_result = readRDS("zeisel_imputation_result.rds")

# Calculate correlation with reference data, CMD between the Pearson correlation 
# matrix derived from the imputed (and observed) matrix and the Pearson correlation 
# matrix derived from the reference matrix
K = length(zeisel_imputation_result$Methods.used)
ncell = dim(zeisel$count.samp)[2]
ngene = dim(zeisel$count.samp)[1]
zeisel_cor.gene = matrix(0,ngene,K+2)
zeisel_cor.cell = matrix(0,ncell,K+2)
zeisel_CMD.gene = rep(0,K+2)
zeisel_CMD.cell = rep(0,K+2)

ref.count =  log(normalizeData(zeisel$count.ref)+1)
ref.cor.gene = cor(t(ref.count))
ref.cor.cell = cor(ref.count)
for (k in 1:K){
  cat("k=",k)
  cat("\n")
  imputed.count = log(normalizeData(zeisel_imputation_result$count.imputed.individual[,,k])+1)
  zeisel_cor.gene[,k] = t(get.cor.gene(ref.count, imputed.count))
  zeisel_cor.cell[,k] = t(get.cor.cell(ref.count, imputed.count))
  imputed.cor.gene = cor(t(imputed.count))
  imputed.cor.gene[is.na(imputed.cor.gene)] = 0
  imputed.cor.cell = cor(imputed.count)
  imputed.cor.cell[is.na(imputed.cor.cell)] = 0
  zeisel_CMD.gene[k] = calc_cmd(ref.cor.gene, imputed.cor.gene)
  zeisel_CMD.cell[k] = calc_cmd(ref.cor.cell, imputed.cor.cell)
}

EnImpute.count =  log(normalizeData(zeisel_imputation_result$count.EnImpute.exp)+1)
zeisel_cor.gene[,K+1] = t(get.cor.gene(ref.count, EnImpute.count))
zeisel_cor.cell[,K+1] = t(get.cor.cell(ref.count, EnImpute.count))
EnImpute.cor.gene = cor(t(EnImpute.count))
EnImpute.cor.cell = cor(EnImpute.count)
zeisel_CMD.gene[K+1] = calc_cmd(ref.cor.gene, EnImpute.cor.gene)
zeisel_CMD.cell[K+1] = calc_cmd(ref.cor.cell, EnImpute.cor.cell)


sample.count =  log(normalizeData(zeisel$count.samp)+1)
zeisel_cor.gene[,K+2] = t(get.cor.gene(ref.count, sample.count))
zeisel_cor.cell[,K+2] = t(get.cor.cell(ref.count, sample.count))
sample.cor.gene = cor(t(sample.count))
sample.cor.cell = cor(sample.count)
zeisel_CMD.gene[K+2] = calc_cmd(ref.cor.gene, sample.cor.gene)
zeisel_CMD.cell[K+2] = calc_cmd(ref.cor.cell, sample.cor.cell)

save(zeisel, zeisel_imputation_result, zeisel_cor.gene, zeisel_cor.cell, zeisel_CMD.gene, zeisel_CMD.cell, file = "zeisel_result.RData")




##################################################
# run clustering and tsne on the baron dataset
load("baron_result.RData")
baron = readRDS("baron.rds")
baron.cluster.ref = tsne_cluster(baron$count.ref, res =1, pcs=15)
baron.cluster.observed = tsne_cluster(baron$count.samp, res =1, pcs=15)
baron.cluster.EnImpute = tsne_cluster(exp(baron_imputation_result$count.EnImpute.log)-1,  res =1, pcs=15)

baron.cluster.imputed = list()
K = length(baron_imputation_result$Methods.used)
for (k in 1:K){
  cat("k=",k)
  baron.cluster.imputed[[k]] =  tsne_cluster(baron_imputation_result$count.imputed.individual[,,k], res =1, pcs=15)
}
save(baron.cluster.imputed, baron.cluster.observed, baron.cluster.ref, baron.cluster.EnImpute, file = "baron.cluster.Rdata")


load("baron.cluster.Rdata")
K = length(baron.cluster.imputed)
# plot text
par(mfrow = c(3, 4),  mar=c(1.5,1.5,1.5,1.5)) 
for (k in 1:K){
  plot(baron.cluster.imputed[[k]]@dr$tsne@cell.embeddings, col = as.numeric(baron.cluster.ref@ident), pch = 19, axes = FALSE, 
       frame.plot = TRUE, ann = FALSE, cex = 0.6)
}
plot(baron.cluster.EnImpute@dr$tsne@cell.embeddings, col = as.numeric(baron.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
plot(baron.cluster.observed@dr$tsne@cell.embeddings, col = as.numeric(baron.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
plot(baron.cluster.ref@dr$tsne@cell.embeddings, col = as.numeric(baron.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)

s1 = 0.12
d1 = 0.25
s2 = -1
d2 = -13.75
mtext(baron_imputation_result$Methods.used[[1]], outer = TRUE, at = s1 , line = s2, cex = 1)
mtext(baron_imputation_result$Methods.used[[2]], outer = TRUE, at = s1  + d1, line = s2, cex = 1)
mtext(baron_imputation_result$Methods.used[[3]], outer = TRUE, at = s1  + 2*d1, line = s2, cex = 1)
mtext(baron_imputation_result$Methods.used[[4]], outer = TRUE, at = s1  + 3*d1, line = s2, cex = 1)
mtext(baron_imputation_result$Methods.used[[5]], outer = TRUE, at = s1 , line = s2+d2, cex = 1)
mtext(baron_imputation_result$Methods.used[[6]], outer = TRUE, at = s1  + d1, line = s2+d2, cex = 1)
mtext(baron_imputation_result$Methods.used[[7]], outer = TRUE, at = s1  + 2*d1, line = s2+d2, cex = 1)
mtext(baron_imputation_result$Methods.used[[8]], outer = TRUE, at = s1  + 3*d1, line = s2+d2, cex = 1)
mtext("EnImpute", outer = TRUE, at = s1, line = s2+ 2*d2, cex = 1)
mtext("Observed", outer = TRUE,  at = s1  + d1, line = s2+ 2*d2, cex =1)
mtext("Reference", outer = TRUE, at = s1  + 2*d1, line = s2+ 2*d2, cex = 1)









##################################################
#  run clustering and tsne on the manno dataset
load("manno_result.RData")
manno = readRDS("manno.rds")
load("manno_result.RData")
manno = readRDS("manno.rds")
manno.cluster.ref = tsne_cluster(manno$count.ref, res =1, pcs=15)
manno.cluster.observed = tsne_cluster(manno$count.samp, res =1, pcs=15)
manno.cluster.EnImpute = tsne_cluster(exp(manno_imputation_result$count.EnImpute.log)-1,  res =1, pcs=15)

manno.cluster.imputed = list()
K = length(manno_imputation_result$Methods.used)

for (k in 1:K){
  cat("k=",k)
  manno.cluster.imputed[[k]] =  tsne_cluster(manno_imputation_result$count.imputed.individual[,,k], res =1, pcs=15)
}
save(manno.cluster.imputed, manno.cluster.observed, manno.cluster.ref, manno.cluster.EnImpute, file = "manno.cluster.Rdata")

load("manno.cluster.Rdata")
K = length(manno.cluster.imputed)
# plot text
par(mfrow = c(3, 4),  mar=c(1.5,1.5,1.5,1.5)) 
for (k in 1:K){
  plot(manno.cluster.imputed[[k]]@dr$tsne@cell.embeddings, col = as.numeric(manno.cluster.ref@ident), pch = 19, axes = FALSE, 
       frame.plot = TRUE, ann = FALSE, cex = 0.6)
}
plot(manno.cluster.EnImpute@dr$tsne@cell.embeddings, col = as.numeric(manno.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
plot(manno.cluster.observed@dr$tsne@cell.embeddings, col = as.numeric(manno.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
plot(manno.cluster.ref@dr$tsne@cell.embeddings, col = as.numeric(manno.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)

s1 = 0.12
d1 = 0.25
s2 = -1
d2 = -13.75
mtext(manno_imputation_result$Methods.used[[1]], outer = TRUE, at = s1 , line = s2, cex = 1)
mtext(manno_imputation_result$Methods.used[[2]], outer = TRUE, at = s1  + d1, line = s2, cex = 1)
mtext(manno_imputation_result$Methods.used[[3]], outer = TRUE, at = s1  + 2*d1, line = s2, cex = 1)
mtext(manno_imputation_result$Methods.used[[4]], outer = TRUE, at = s1  + 3*d1, line = s2, cex = 1)
mtext(manno_imputation_result$Methods.used[[5]], outer = TRUE, at = s1 , line = s2+d2, cex = 1)
mtext(manno_imputation_result$Methods.used[[6]], outer = TRUE, at = s1  + d1, line = s2+d2, cex = 1)
mtext(manno_imputation_result$Methods.used[[7]], outer = TRUE, at = s1  + 2*d1, line = s2+d2, cex = 1)
mtext(manno_imputation_result$Methods.used[[8]], outer = TRUE, at = s1  + 3*d1, line = s2+d2, cex = 1)
mtext("EnImpute", outer = TRUE, at = s1, line = s2+ 2*d2, cex = 1)
mtext("Observed", outer = TRUE,  at = s1  + d1, line = s2+ 2*d2, cex =1)
mtext("Reference", outer = TRUE, at = s1  + 2*d1, line = s2+ 2*d2, cex = 1)





##################################################
# run clustering and tsne on the chen dataset
load("chen_result.RData")
chen = readRDS("chen.rds")
load("chen_result.RData")
chen = readRDS("chen.rds")
chen.cluster.ref = tsne_cluster(chen$count.ref, res =1, pcs=15)
chen.cluster.observed = tsne_cluster(chen$count.samp, res =1, pcs=15)
chen.cluster.EnImpute = tsne_cluster(exp(chen_imputation_result$count.EnImpute.log)-1,  res =1, pcs=15)

chen.cluster.imputed = list()
K = length(chen_imputation_result$Methods.used)

for (k in 1:K){
  cat("k=",k)
  chen.cluster.imputed[[k]] =  tsne_cluster(chen_imputation_result$count.imputed.individual[,,k], res =1, pcs=15)
}
save(chen.cluster.imputed, chen.cluster.observed, chen.cluster.ref, chen.cluster.EnImpute, file = "chen.cluster.Rdata")

load("chen.cluster.Rdata")
K = length(chen.cluster.imputed)
# plot text
par(mfrow = c(3, 4),  mar=c(1.5,1.5,1.5,1.5)) 
for (k in 1:K){
  plot(chen.cluster.imputed[[k]]@dr$tsne@cell.embeddings, col = as.numeric(chen.cluster.ref@ident), pch = 19, axes = FALSE, 
       frame.plot = TRUE, ann = FALSE, cex = 0.6)
}
plot(chen.cluster.EnImpute@dr$tsne@cell.embeddings, col = as.numeric(chen.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
plot(chen.cluster.observed@dr$tsne@cell.embeddings, col = as.numeric(chen.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
plot(chen.cluster.ref@dr$tsne@cell.embeddings, col = as.numeric(chen.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)

s1 = 0.12
d1 = 0.25
s2 = -1
d2 = -13.75
mtext(chen_imputation_result$Methods.used[[1]], outer = TRUE, at = s1 , line = s2, cex = 1)
mtext(chen_imputation_result$Methods.used[[2]], outer = TRUE, at = s1  + d1, line = s2, cex = 1)
mtext(chen_imputation_result$Methods.used[[3]], outer = TRUE, at = s1  + 2*d1, line = s2, cex = 1)
mtext(chen_imputation_result$Methods.used[[4]], outer = TRUE, at = s1  + 3*d1, line = s2, cex = 1)
mtext(chen_imputation_result$Methods.used[[5]], outer = TRUE, at = s1 , line = s2+d2, cex = 1)
mtext(chen_imputation_result$Methods.used[[6]], outer = TRUE, at = s1  + d1, line = s2+d2, cex = 1)
mtext(chen_imputation_result$Methods.used[[7]], outer = TRUE, at = s1  + 2*d1, line = s2+d2, cex = 1)
mtext("EnImpute", outer = TRUE, at = s1 + 3*d1, line = s2+ d2, cex = 1)
mtext("Observed", outer = TRUE,  at = s1, line = s2+ 2*d2, cex =1)
mtext("Reference", outer = TRUE, at = s1  + d1, line = s2+ 2*d2, cex = 1)





##################################################
# run clustering and tsne on the zeisel dataset
load("zeisel_result.RData")
zeisel = readRDS("zeisel.rds")
load("zeisel_result.RData")
zeisel = readRDS("zeisel.rds")
zeisel.cluster.ref = tsne_cluster(zeisel$count.ref, res =1, pcs=15)
zeisel.cluster.observed = tsne_cluster(zeisel$count.samp, res =1, pcs=15)
zeisel.cluster.EnImpute = tsne_cluster(exp(zeisel_imputation_result$count.EnImpute.log)-1,  res =1, pcs=15)


zeisel.cluster.imputed = list()
K = length(zeisel_imputation_result$Methods.used)

for (k in 1:K){
  cat("k=",k)
  zeisel.cluster.imputed[[k]] =  tsne_cluster(zeisel_imputation_result$count.imputed.individual[,,k], res =1, pcs=15)
}
save(zeisel.cluster.imputed, zeisel.cluster.observed, zeisel.cluster.ref, zeisel.cluster.EnImpute, file = "zeisel.cluster.Rdata")



load("zeisel.cluster.Rdata")
K = length(zeisel.cluster.imputed)
# plot text
par(mfrow = c(3, 4),  mar=c(1.5,1.5,1.5,1.5)) 
for (k in 1:K){
  plot(zeisel.cluster.imputed[[k]]@dr$tsne@cell.embeddings, col = as.numeric(zeisel.cluster.ref@ident), pch = 19, axes = FALSE, 
       frame.plot = TRUE, ann = FALSE, cex = 0.6)
}
plot(zeisel.cluster.EnImpute@dr$tsne@cell.embeddings, col = as.numeric(zeisel.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
plot(zeisel.cluster.observed@dr$tsne@cell.embeddings, col = as.numeric(zeisel.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
plot(zeisel.cluster.ref@dr$tsne@cell.embeddings, col = as.numeric(zeisel.cluster.ref@ident), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)

s1 = 0.12
d1 = 0.25
s2 = -1
d2 = -13.75
mtext(zeisel_imputation_result$Methods.used[[1]], outer = TRUE, at = s1 , line = s2, cex = 1)
mtext(zeisel_imputation_result$Methods.used[[2]], outer = TRUE, at = s1  + d1, line = s2, cex = 1)
mtext(zeisel_imputation_result$Methods.used[[3]], outer = TRUE, at = s1  + 2*d1, line = s2, cex = 1)
mtext(zeisel_imputation_result$Methods.used[[4]], outer = TRUE, at = s1  + 3*d1, line = s2, cex = 1)
mtext(zeisel_imputation_result$Methods.used[[5]], outer = TRUE, at = s1 , line = s2+d2, cex = 1)
mtext(zeisel_imputation_result$Methods.used[[6]], outer = TRUE, at = s1  + d1, line = s2+d2, cex = 1)
mtext(zeisel_imputation_result$Methods.used[[7]], outer = TRUE, at = s1  + 2*d1, line = s2+d2, cex = 1)
mtext(zeisel_imputation_result$Methods.used[[8]], outer = TRUE, at = s1  + 3*d1, line = s2+d2, cex = 1)
mtext("EnImpute", outer = TRUE, at = s1, line = s2+ 2*d2, cex = 1)
mtext("Observed", outer = TRUE,  at = s1  + d1, line = s2+ 2*d2, cex =1)
mtext("Reference", outer = TRUE, at = s1  + 2*d1, line = s2+ 2*d2, cex = 1)





###################################################################3
## Compare clustering with reference
rm(list=ls())
library(mclust)
sim = list();
load("baron.cluster.Rdata")
baron_imputation_result = readRDS("baron_imputation_result.rds")
sim[[1]] = rep(NA,10);
K =length(baron_imputation_result$Methods.used)
for (k in 1:K){
  sim[[1]][k] = adjustedRandIndex(baron.cluster.imputed[[k]]@ident, baron.cluster.ref@ident)
}
sim[[1]][K+1] = adjustedRandIndex(baron.cluster.ref@ident, baron.cluster.EnImpute@ident)
sim[[1]][K+2] = adjustedRandIndex(baron.cluster.ref@ident, baron.cluster.observed@ident)
names(sim[[1]]) = c(baron_imputation_result$Methods.used,"EnImpute","observed")
rm(baron.cluster.EnImpute, baron.cluster.imputed, baron.cluster.observed, baron.cluster.ref, baron_imputation_result)

load("chen.cluster.Rdata")
chen_imputation_result = readRDS("chen_imputation_result.rds")
sim[[2]] = rep(NA,9);
K =length(chen_imputation_result$Methods.used)
for (k in 1:K){
  sim[[2]][k] = adjustedRandIndex(chen.cluster.imputed[[k]]@ident, chen.cluster.ref@ident)
}
sim[[2]][K+1] = adjustedRandIndex(chen.cluster.ref@ident, chen.cluster.EnImpute@ident)
sim[[2]][K+2] = adjustedRandIndex(chen.cluster.ref@ident, chen.cluster.observed@ident)
names(sim[[2]]) = c(chen_imputation_result$Methods.used,"EnImpute","observed")
rm(chen.cluster.EnImpute, chen.cluster.imputed, chen.cluster.observed, chen.cluster.ref, chen_imputation_result)

load("manno.cluster.Rdata")
manno_imputation_result = readRDS("manno_imputation_result.rds")
sim[[3]] = rep(NA,10);
K =length(manno_imputation_result$Methods.used)
for (k in 1:K){
  sim[[3]][k] = adjustedRandIndex(manno.cluster.imputed[[k]]@ident, manno.cluster.ref@ident)
}
sim[[3]][K+1] = adjustedRandIndex(manno.cluster.ref@ident, manno.cluster.EnImpute@ident)
sim[[3]][K+2] = adjustedRandIndex(manno.cluster.ref@ident, manno.cluster.observed@ident)
names(sim[[3]]) = c(manno_imputation_result$Methods.used,"EnImpute","observed")
rm(manno.cluster.EnImpute, manno.cluster.imputed, manno.cluster.observed, manno.cluster.ref, manno_imputation_result)


load("zeisel.cluster.Rdata")
zeisel_imputation_result = readRDS("zeisel_imputation_result.rds")
sim[[4]] = rep(NA,10);
K =length(zeisel_imputation_result$Methods.used)
for (k in 1:K){
  sim[[4]][k] = adjustedRandIndex(zeisel.cluster.imputed[[k]]@ident, zeisel.cluster.ref@ident)
}
sim[[4]][K+1] = adjustedRandIndex(zeisel.cluster.ref@ident, zeisel.cluster.EnImpute@ident)
sim[[4]][K+2] = adjustedRandIndex(zeisel.cluster.ref@ident, zeisel.cluster.observed@ident)
names(sim[[4]]) = c(zeisel_imputation_result$Methods.used,"EnImpute","observed")
rm(zeisel.cluster.EnImpute, zeisel.cluster.imputed, zeisel.cluster.observed, zeisel.cluster.ref, zeisel_imputation_result)


ll = 0.3
ee = 0.1
par(mfrow = c(4, 1),  mar=c(ee,ee,ee,ee), mai=c(ll,ll,ll,ll)) 
barplot(sim[[1]], col = rainbow(10))
tmp =rep(NA,10)
tmp[1:2] = sim[[2]][1:2]
tmp[4:10] = sim[[2]][3:9]
names(tmp) = names(sim[[1]])
barplot(tmp, col = rainbow(10))
barplot(sim[[3]], col = rainbow(10))
barplot(sim[[4]], col = rainbow(10))




###################################################################3
## Correlation with reference (cell)
load("baron_result.RData")
load("chen_result.RData")
load("manno_result.RData")
load("zeisel_result.RData")

ll = 0.3
ee = 0.1
par(mfrow = c(4, 1),  mar=c(ee,ee,ee,ee), mai=c(ll,ll,ll,ll)) 
label = c(baron_imputation_result$Methods.used,"EnImpute", "Observed")
colnames(baron_cor.cell) = label
boxplot(baron_cor.cell, col = rainbow(10),  notch=T, ylim=c(0.3,0.7))
tmp = matrix(NA, dim(chen_cor.cell)[1],  dim(chen_cor.cell)[2]+1)
tmp[,1:2] = chen_cor.cell[,1:2]
tmp[,4:10] = chen_cor.cell[,3:9]
colnames(tmp) = label
boxplot(tmp, col = rainbow(10),  notch=T, ylim=c(0.2,0.7))
colnames(manno_cor.cell) = label
boxplot(manno_cor.cell, col = rainbow(10), notch=TRUE, ylim=c(0.3,0.7))
colnames(zeisel_cor.cell) = label
boxplot(zeisel_cor.cell, col = rainbow(10), notch=TRUE, ylim=c(0.3,0.7))



###################################################################3
## Correlation with reference (gene)
load("baron_result.RData")
load("chen_result.RData")
load("manno_result.RData")
load("zeisel_result.RData")

ll = 0.3
ee = 0.1
par(mfrow = c(4, 1),  mar=c(ee,ee,ee,ee), mai=c(ll,ll,ll,ll)) 
label = c(baron_imputation_result$Methods.used,"EnImpute", "Observed")
colnames(baron_cor.gene) = label
boxplot(baron_cor.gene, col = rainbow(10), notch=TRUE, ylim = c(0.1,0.5))
tmp = matrix(NA, dim(chen_cor.gene)[1],  dim(chen_cor.gene)[2]+1)
tmp[,1:2] = chen_cor.gene[,1:2]
tmp[,4:10] = chen_cor.gene[,3:9]
colnames(tmp) = label
boxplot(tmp, col = rainbow(10),  notch=T, ylim = c(0.1,0.5))
colnames(manno_cor.gene) = label
boxplot(manno_cor.gene, col = rainbow(10), notch=TRUE, ylim = c(0.1,0.5))
colnames(zeisel_cor.gene) = label
boxplot(zeisel_cor.gene, col = rainbow(10), notch=TRUE, ylim = c(0.1,0.6))




###################################################################
## Comparison of  cell-to-cell correlation matrices of recovered values 
# with the true correlation matrices
load("baron_result.RData")
load("chen_result.RData")
load("manno_result.RData")
load("zeisel_result.RData")
ll = 0.3
ee = 0.1
par(mfrow = c(4, 1),  mar=c(ee,ee,ee,ee), mai=c(ll,ll,ll,ll)) 

label = c(baron_imputation_result$Methods.used,"EnImpute", "Observed")
names(baron_CMD.cell) = label
barplot(baron_CMD.cell, col = rainbow(10))
tmp =rep(NA,10)
tmp[1:2] = chen_CMD.cell[1:2]
tmp[4:10] = chen_CMD.cell[3:9]
names(tmp) = label
barplot(tmp, col = rainbow(10))
names(manno_CMD.cell) = label
barplot(manno_CMD.cell, col = rainbow(10))
names(zeisel_CMD.cell) = label
barplot(zeisel_CMD.cell, col = rainbow(10))



###################################################################
## Comparison of  gene-to-gene correlation matrices of recovered values 
# with the true correlation matrices
load("baron_result.RData")
load("chen_result.RData")
load("manno_result.RData")
load("zeisel_result.RData")
ll = 0.3
ee = 0.1
par(mfrow = c(4, 1),  mar=c(ee,ee,ee,ee), mai=c(ll,ll,ll,ll)) 

label = c(baron_imputation_result$Methods.used,"EnImpute", "Observed")
names(baron_CMD.gene) = label
barplot(baron_CMD.gene, col = rainbow(10))
tmp =rep(NA,10)
tmp[1:2] = chen_CMD.gene[1:2]
tmp[4:10] = chen_CMD.gene[3:9]
names(tmp) = label
barplot(tmp, col = rainbow(10))
names(manno_CMD.gene) = label
barplot(manno_CMD.gene, col = rainbow(10))
names(zeisel_CMD.gene) = label
barplot(zeisel_CMD.gene, col = rainbow(10))

