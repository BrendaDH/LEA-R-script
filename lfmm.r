library(LEA)
### sNMF
project <- NULL
project <- snmf("file.geno", K=1:15, entropy = TRUE, repetitions = 10, project = "new", iterations = 100000)
# plot cross-entropy criterion of all runs of the project
plot(project, cex = 1.2, col = "lightblue", pch = 19)
# get the cross-entropy of each run for K = 8 (e.g. 169 samples, all species)
ce <- cross.entropy(project, K = 8)
# select the run with the lowest cross-entropy
best <- which.min(ce)
best

### LFMM
obj.lfmm <- lfmm("file.lfmm", "file.env", K = 8, rep = 5, iterations = 200000, burnin = 50000, project="new")

# correcting p-values
qv1 = which(qvalue(adj.p.values1, fdr = .1)$signif)
zs1 = z.scores(obj.lfmm, K = 8, d=1)
zs.median1 = apply(zs1, MARGIN = 1, median)
adj.p.values1 = pchisq(zs.median1^2/1, df = 1, lower = FALSE)
qv1 = which(qvalue(adj.p.values1, fdr = .1)$signif)

par(mfrow=c(3,3))
# histogram
hist(adj.p.values1, col = "orangered")
# Manhattan plot
plot(-log10(adj_pvals1), pch = 19, col = "royalblue1", cex = .7, ylab="-log10(pvalues)")

# and so on replacing argument d= from 2 to 13