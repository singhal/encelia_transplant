library(ggplot2)
library(cowplot)
library(LEA)
library(adegenet)

a = read.table("~/Dropbox/Encelia/analysis/hybrid_zones/admixture/ASPVEN.miss0.5.2.Q")
names(a) = c("admix_asp", "admix_ven")
n = read.table("~/Dropbox/Encelia/analysis/hybrid_zones/NGSadmix/ASPVEN.2.qopt")
names(n) = c("ngs_ven", "ngs_asp")
an = cbind(a, n)
figB = ggplot(an, aes(admix_asp, ngs_asp)) + 
  geom_point() + xlab("ADMIXTURE, E. asperifolia ancestry") +
  ylab("NGSadmix, E. asperifolia ancestry") +
  geom_abline(intercept = 0, slope = 1, col = "Red")
cor.test(an$admix_asp, an$ngs_asp)


a = read.table("~/Dropbox/Encelia/analysis/hybrid_zones/admixture/PALVEN.miss0.5.2.Q")
names(a) = c("admix_ven", "admix_pal")
n = read.table("~/Dropbox/Encelia/analysis/hybrid_zones/NGSadmix/PALVEN.2.qopt")
names(n) = c("ngs_pal", "ngs_ven")
an = cbind(a, n)
figA = ggplot(an, aes(admix_pal, ngs_pal)) + 
  geom_point() + xlab("ADMIXTURE, E. palmeri ancestry") +
  ylab("NGSadmix, E. palmeri ancestry") +
  geom_abline(intercept = 0, slope = 1, col = "Red")
cor.test(an$admix_pal, an$ngs_pal)


ab = (figA | figB) +  plot_annotation(tag_levels = 'A')
save_plot("~/Dropbox/Encelia/manuscripts/EnceliaTransplant_v2/figures/ancestry_estimates_across_programs.png", 
          ab, nrow = 2, base_width = 8, base_height = 2)
