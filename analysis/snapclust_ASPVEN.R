library(adegenet)
library(dartR)
library(cowplot)
theme_set(theme_cowplot())

hybtypes = c("A", "B", "F1", "F2", "bcA1", "bcB1",
             "bcA2", "bcB2", "bcA3", "bcB3")
names(hybtypes) = c("A", "B", "0.5_A-0.5_B", "F2", "0.75_A-0.25_B",
                    "0.25_A-0.75_B", "0.875_A-0.125_B", 
                    "0.125_A-0.875_B", "0.0163_B-0.9837_A",
                    "0.0163_A-0.9837_B")

# load in SNP data
d = read.snp("~/Dropbox/Encelia/analysis/hybrid_zones/input_files/ASPVEN.miss0.8.snp")

# d1 = gl.subsample.loci(d, 1000, method = "random")
d1 = as.matrix(d)
d1[which(d1 == 0)] <- "1/1" # homozygote reference
d1[which(d1 == 1)] <- "1/2" # heterozygote
d1[which(d1 == 2)] <- "2/2" # homozygote alternate
d2 <- df2genind(d1, sep = "/", ploidy = 2)

grp.ini <- find.clusters(d2, n.clust=2, n.pca=15)
resh <- snapclust(d2, 2, pop.ini = grp.ini$grp,
                  hybrids = TRUE, hybrid.coef = c(0.01625, 0.125, .25, .5))
probs = resh$proba
res = data.frame(ind = rownames(resh$proba),
                 hybtype = colnames(probs)[ apply(probs, 1, function(x) { which(x == max(x)) }) ],
                 stringsAsFactors = F)
res$hybtype = hybtypes[res$hybtype]

tbl = read.table("~/Dropbox/Encelia/analysis/hybrid_zones/NGSadmix/ASPVEN.2.qopt",
                 col.names = c("Ancestry1", "Ancestry2"))
inds = read.csv("~/Dropbox/Encelia/analysis/hybrid_zones/ASPVEN_inds.txt", header=F, stringsAsFactors = F)
tbl$sample = inds$V1

res = inner_join(res, tbl, by = c("ind" = "sample"))
saveRDS(res, "~/Desktop/ASPVEN_hybclass")

#######################
# do hybridization
#######################

p1ix = which(colnames(probs)[ apply(probs, 1, function(x) { which(x == max(x)) }) ] == "A")
p2ix = which(colnames(probs)[ apply(probs, 1, function(x) { which(x == max(x)) }) ] == "B")

d2a = d2[p1ix, ]
pop(d2a) = rep("A", nInd(d2a))
d2b = d2[p2ix, ]
pop(d2b) = rep("B", nInd(d2b))

f1 <- hybridize(d2a, d2b, n = 50, pop = "F1")
indNames(f1) = sapply(seq(1:nInd(f1)), function(x) {paste0("f1_", x)})

# back cross gen 1 - 0.25
bcA1 <- hybridize(f1, d2a, n = 50, pop = "bcA1")
indNames(bcA1) = sapply(seq(1:nInd(bcA1)), function(x) {paste0("bcA1_", x)})
bcB1 <- hybridize(f1, d2b, n = 50, pop = "bcB1")
indNames(bcB1) = sapply(seq(1:nInd(bcB1)), function(x) {paste0("bcB1_", x)})

# back cross gen 2 - 0.125
bcA2 <- hybridize(bcA1, d2a, n = 50, pop = "bcA2")
indNames(bcA2) = sapply(seq(1:nInd(bcA2)), function(x) {paste0("bcA2_", x)})
bcB2 <- hybridize(bcB1, d2b, n = 50, pop = "bcB2")
indNames(bcB2) = sapply(seq(1:nInd(bcB2)), function(x) {paste0("bcB2_", x)})

# back cross gen 2 - 0.0625
bcA3 <- hybridize(bcA2, d2a, n = 50, pop = "bcA3")
indNames(bcA3) = sapply(seq(1:nInd(bcA3)), function(x) {paste0("bcA3_", x)})
bcB3 <- hybridize(bcB2, d2b, n = 50, pop = "bcB3")
indNames(bcB3) = sapply(seq(1:nInd(bcB3)), function(x) {paste0("bcB3_", x)})


f2 <- hybridize(f1, f1, n = 50, pop = "F2")
indNames(f2) = sapply(seq(1:nInd(f2)), function(x) {paste0("f2_", x)})
x <- repool(d2a, d2b, f1, f2, bcA1, bcA2, bcA3, bcB1, bcB2, bcB3)
ressim <- snapclust(x, 2, hybrids = TRUE, 
                    hybrid.coef = c(0.01625, 0.125, .25, .5))

infer = colnames(ressim$proba)[ apply(ressim$proba, 1, function(x) { which(x == max(x)) }) ]

df <- data.frame(table(pop(x), infer))
df$infer2 = hybtypes[as.character(df$infer)]

df$Var1 <- factor(df$Var1, levels = c(hybtypes))
df$infer2 <- factor(df$infer2, levels = c(hybtypes))

saveRDS(df, "~/Desktop/ASPVEN_sim_hybclass")
