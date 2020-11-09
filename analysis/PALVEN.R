library(ape)
library(geosphere)
library(hzar)
library(readr)
library(cowplot)
library(ggplot2)
library(sp)
library(stringr)
library(dplyr)
theme_set(theme_cowplot())

###########################
# PALMERI - VENTORUM HZ
###########################

inds1 = read_csv("~/Dropbox/Encelia/leaf_images/organized/image_to_sample_map.csv")
inds2 = read_csv("~/Dropbox/Encelia/analysis/hybrid_zones/Encelia_Samples - GENERAL.csv")
# leaf data
p1 = read_csv("~/Dropbox/Encelia/analysis/hybrid_zones/leaf_images/leaf_image - john.csv")
p1 = p1 %>% dplyr::filter(m != 'IMG_1008')
p2 = left_join(p1, inds1, by = c("m" = "PICTURE_NAME"))

# a few inds are doubled, so merge by sample name not picture name
vals = c("Area", "Mean", "Mode", "Perimeter",
         "Major", "Minor", "Angle", "Circ.", 
         "AR", "Round", "Solidty")
p3 = p2 %>% group_by(SAMPLE_ID, LeafID) %>% 
  dplyr::select(vals) %>% summarize_all(sum, na.rm = T) %>% ungroup()
write.csv(p3, "~/Desktop/PALVEN_leaf_image_data.csv", row.names = F,
          quote = F)

# do the pca
pca2 = prcomp(p3 %>% dplyr::select(vals), center = T, scale. = T)
p4 = cbind(p3, pca2$x)
p5 = p4 %>% group_by(SAMPLE_ID) %>% 
  dplyr::select("PC1", "PC2") %>% summarize_all(mean, na.rm = T) %>% ungroup()

# read in data & match up to ind names
tbl = read.table("~/Dropbox/Encelia/analysis/hybrid_zones/NGSadmix/PALVEN.2.qopt",
                 col.names = c("Ancestry1", "Ancestry2"))
inds = read.csv("~/Dropbox/Encelia/analysis/hybrid_zones/PALVEN_inds.txt", header=F, stringsAsFactors = F)
tbl$sample = inds$V1
tbl = left_join(tbl, p5, by = c("sample" = "SAMPLE_ID"))

# add in lats and longs
tbl$sample2 = gsub("-","", tbl$sample)
tbl = left_join(tbl, 
          inds2 %>% dplyr::select(PLANT_ID, PUTATIVE_SPECIES, LATITUDE, LONGITUDE), 
          by = c("sample2" = "PLANT_ID"))

# add in habitat
hab = read.csv("~/Dropbox/scripts/encelia/figures/palven_habitat.csv", stringsAsFactors = F)
tbl$habitat = hab[match(tbl$sample, hab$SAMPLE_ID), "habitat"]
tbl$habitat = factor(tbl$habitat, ordered = TRUE, 
                     levels = c("desert", "ecotone", "dune"))

# create a posthoc transect
# use E ventorum because mostly linear
ev = tbl[tbl$PUTATIVE_SPECIES == 'Encelia_ventorum', ]
trans = lm(ev$LATITUDE ~ ev$LONGITUDE)

# diagnostic plot
cols = brewer.pal(3, "Set1")
transplot = ggplot(tbl, aes(LONGITUDE, LATITUDE, fill = habitat)) +
  geom_point(shape = 21) + scale_fill_manual(values = cols) +
  geom_abline(slope = trans$coefficients[2], intercept = trans$coefficients[1])
save_plot("~/Desktop/palmeri_ventorum_trans.pdf", transplot, 
          base_height = 5, base_width = 6)

# get distance
lon = seq(min(tbl$LONGITUDE), max(tbl$LONGITUDE), 0.00001)
lat = lon * coef(trans)[2] + coef(trans)[1]
tpoints = dist2Line(as.matrix(tbl[, c("LONGITUDE", "LATITUDE")]), 
                    as.matrix(cbind(lon, lat)))

# find starting point
start = tpoints[which(tpoints[, 3] == max(tpoints[ ,3])), c("lon", "lat")]
dists = rep(NA, nrow(tpoints))
for (i in 1:length(dists)) {
  dists[i] = distHaversine(start, tpoints[i, 2:3])
}
tbl$distance = dists
tbl$PUTATIVE_SPECIES = factor(tbl$PUTATIVE_SPECIES, 
                              levels(as.factor(tbl$PUTATIVE_SPECIES))[c(1, 3, 2)])
tbl = tbl[order(tbl$PUTATIVE_SPECIES, -tbl$distance), ]

# transect distance
max(tbl$distance)
# ecotone transect
eco = tbl[tbl$habitat == "ecotone", ]
minpt = eco[which(eco$LONGITUDE == max(eco$LONGITUDE)), c("LONGITUDE", "LATITUDE")][1, ]
maxpt = eco[which(eco$LONGITUDE == min(eco$LONGITUDE)), c("LONGITUDE", "LATITUDE")][1, ]
distHaversine(minpt, maxpt)

# fit dna cline
# remove  missing data
tbl4 = tbl[complete.cases(tbl$Ancestry2), ]

tbl4$distance2 = round(tbl4$distance / 50, 0) * 50
tbl5 = tbl4 %>% group_by(distance2) %>% summarize(ancestry = mean(Ancestry2)) %>% ungroup()
cl2 = hzar.doMolecularData1DPops(tbl5$distance2, tbl5$ancestry, rep(10, nrow(tbl5)))
clModel2 <- hzar.makeCline1DFreq(cl2, scaling="fixed",tails="none");

clModel2 <- hzar.model.addBoxReq(clModel2, 0, 500);
cl_model_fitR2 = hzar.first.fitRequest.old.ML(clModel2, cl2, verbose=TRUE)

cl_model_fitR2$mcmcParam$chainLength <- 1e5;
cl_model_fitR2$mcmcParam$burnin <- 1e2;
clModelFit22 <- hzar.doFit(cl_model_fitR2)
clModelFitData2 <- hzar.dataGroup.add(clModelFit22)
clModelFitData2 <- hzar.dataGroup.add(clModelFitData2, 
                                      hzar.chain.doSeq(hzar.next.fitRequest(clModelFit22)));

pts = seq(min(tbl5$distance2), max(tbl5$distance2), 1)
fzCoor <- hzar.getCredParamRed(clModelFitData2)$fzCline(pts)
predval = clModelFitData2$ML.cline$clineFunc(pts)

##################################
# plot cline
##################################
plot_cline <- function() {
  dd1 = data.frame(x = fzCoor$x, y = fzCoor$yMin)
  dd2 = data.frame(x = rev(fzCoor$x), y = rev(fzCoor$yMax))
  pp1 = rbind(dd1, dd2)
  cl1 = data.frame(x = pts, y = predval)

  ggplot(tbl4, aes(distance, Ancestry2)) + geom_point(alpha = 0.7, color = "gray30") + 
    xlab("transect distance (m)") + ylab(expression(paste("% ", italic("E. ventorum"),  " ancestry"))) +
    geom_polygon(data = pp1, aes(x, y), fill = "gray30", alpha = 0.5) + 
    geom_line(data = cl1, aes(x, y)) + theme_cowplot(font_size = 12)
}

################################
# plot structure
################################
plot_str <- function() {
  test = tidyr::gather(tbl %>% dplyr::select(sample, Ancestry1, Ancestry2), 
                       key = "species", value = "proportion", Ancestry1, Ancestry2)
  test2 = left_join(test, tbl)
  indorder = pull(test2 %>% dplyr::filter(species == "Ancestry2") %>% 
                    arrange(habitat, Ancestry2) %>% dplyr::select(sample))
  test2$sample2 = factor(test2$sample, levels = indorder)
  
  ggplot(test2) + 
    geom_bar(aes(x = sample2, y = proportion, fill = species), width = 1, size = 0,
             stat = "identity") + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(),
    ) +
    scale_fill_manual(values = c("#e41a1c", "#377eb8"))
}

#######################
# plot correlation
#######################
plot_corr <- function() {
    ggplot(tbl, aes(PC1, Ancestry2, fill = Ancestry2)) +
    xlab("leaf shape, PC1") + 
    ylab(expression(paste("% ", italic(" E. ventorum"), " ancestry"))) + 
    geom_point( pch = 21) +
    scale_fill_gradient2(high = "#377eb8", mid = "white", low = "#e41a1c", midpoint = 0.5) + 
    theme_cowplot(font_size = 12) + theme(legend.position = "none")
}


all = plot_grid(plot_cline(), plot_corr(), plot_str(), ncol = 1, nrow = 3, align = "h", axis = "l")
save_plot("~/Desktop/PALVEN.pdf", all, ncol = 1, nrow = 3, base_height = 2.5, base_width = 3)


#######################
# plot map
#######################

xx = ggplot(tbl, aes(LONGITUDE, LATITUDE)) + 
  geom_point(aes(shape = as.factor(habitat), fill = Ancestry1),
             position = position_jitter(w = 0.00002, h = 0.00002)) +
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", midpoint = 0.5,
                       name = "E. palmeri ancestry") +
  scale_shape_manual(values = c(21, 22, 24), name = "habitat")
save_plot("~/Desktop/palven.pdf", xx, base_height = 6, base_width = 8)
