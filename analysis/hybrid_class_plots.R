library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
theme_set(theme_cowplot())

setwd("~/Dropbox/Encelia/analysis/hybrid_zones/hybrid_class/")

hybtypes = c("A", "B", "F1", "F2", "bcA1", "bcB1",
             "bcA2", "bcB2", "bcA3", "bcB3")
frac = c(1, 0, 0.5, 0.5, 0.75, 0.25, 
         0.875, 0.125, 0.9375, 0.0625)
exp = data.frame(hybtypes = factor(hybtypes, levels = hybtypes),
                 Ancestry1 = frac)

############
# boxplots
###########

pv = readRDS("PALVEN_hybclass")
pv = pv[complete.cases(pv$hybtype), ]
av = readRDS("ASPVEN_hybclass")

pv$hybtype = factor(pv$hybtype, levels = hybtypes)
av$hybtype = factor(av$hybtype, levels = hybtypes)

a = ggplot(pv, aes(hybtype, Ancestry1)) + geom_boxplot() + 
  xlab("hybrid type") + ylab("ancestry of A") +
  geom_point(data = exp, aes(hybtypes, Ancestry1), col = "red", 
             shape = 8)
b = ggplot(av, aes(hybtype, Ancestry1)) + geom_boxplot() + 
  xlab("hybrid type") + ylab("ancestry of A") +
  geom_point(data = exp, aes(hybtypes, Ancestry1), col = "red", 
             shape = 8)

ab = (a / b) +  plot_annotation(tag_levels = 'A')
save_plot("~/Desktop/ancestry_by_hybridclass.png", ab, nrow = 2)

############
# hyb habitat
###########

# get habitat data for pv
hab = read.csv("~/Dropbox/scripts/encelia/figures/palven_habitat.csv", stringsAsFactors = F)
hab$SAMPLE_ID = gsub("-", "", hab$SAMPLE_ID)
pv2 = left_join(pv, hab, by = c("ind" = "SAMPLE_ID")) %>% 
  dplyr::select(ind, hybtype, habitat, Ancestry1, Ancestry2)
pv2$hybzone = "PALVEN"

# get habitat data for av
tbl = read.csv("~/Dropbox/Encelia/analysis/hybrid_zones/Encelia_Samples - GENERAL.csv", stringsAsFactors = F)
tbl$habitat = NA
tbl[grep("desert", tbl$NOTES), "habitat"] = "desert"
tbl[grep("at ecotone", tbl$NOTES), "habitat"] = "ecotone"
tbl[grep("further", tbl$NOTES), "habitat"] = "dune"
tbl[grep("arroyo", tbl$NOTES), "habitat"] = "ecotone"
av2 = left_join(av, tbl, by = c("ind" = "PLANT_ID"))

# these aren't part of the hybrid zone
outliers = c("ASPVEN11", "ASPVEN111", 
             "ASPVEN112", "ASPVEN113", 
             "ASPVEN114", "ASPVEN115",
             "ASPVEN116", "ASPVEN117",
             "ASPVEN118", "ASPVEN119", 
             "ASPVEN120", "ASPVEN121")
av2 = av2 %>% filter(!ind %in% outliers) %>%
  dplyr::select(ind, hybtype, habitat,  Ancestry1, Ancestry2)
av2$hybzone = "ASPVEN"

ll3 = rbind(pv2, av2)
ll3$hybtype = as.character(ll3$hybtype)

ll3[which(ll3$hybtype == "A" & ll3$hybzone == "ASPVEN"), "hybtype"] = "ventorum"
ll3[which(ll3$hybtype == "B" & ll3$hybzone == "ASPVEN"), "hybtype"] = "asperifolia"
ll3[which(ll3$hybtype == "A" & ll3$hybzone == "PALVEN"), "hybtype"] = "palmeri"
ll3[which(ll3$hybtype == "B" & ll3$hybzone == "PALVEN"), "hybtype"] = "ventorum"
ll3[grep("bc", ll3$hybtype), "hybtype"] = "backcross"


ll3$habitat = factor(ll3$habitat, ordered = TRUE, 
                     levels = c("desert", "ecotone", "dune"))
ll3$species = factor(ll3$hybtype,  
                     levels = c("asperifolia", "palmeri", 
                                "F1", "backcross", "ventorum"))

g1 = ggplot(ll3 %>% filter(complete.cases(habitat))) + 
  geom_bar(aes(x = habitat, fill = species), position = "fill", color = "black") + 
  facet_grid(~hybzone) + 
  scale_fill_manual(values = c("white", "gray90", "gray60", 
                               "gray30", "black"),
                    labels = c(expression(italic("E. asperifolia")),
                               expression(italic("E. palmeri")),
                               "F1", "backcross",
                               expression(italic("E. ventorum")))) +
  theme_cowplot() + theme(legend.text.align = 0)
save_plot("~/Desktop/hyrids_across_habitats.png", g1,
          base_width = 8, base_height = 5)

table(ll3$hybtype, ll3$habitat, ll3$hybzone)

##############
# hybrid index plot
##############

ll3$ventorum = NA
ll3[ll3$hybzone == "PALVEN", "ventorum"] = ll3[ll3$hybzone == "PALVEN", "Ancestry2"]
ll3[ll3$hybzone == "ASPVEN", "ventorum"] = ll3[ll3$hybzone == "ASPVEN", "Ancestry1"]
ggplot(ll3, aes(ventorum)) + 
  geom_histogram(bins = 12) + 
  facet_grid(vars(hybzone), vars(habitat)) + 
  xlab("percent E. ventorum ancestry") + theme_bw()

############
# simulations
###########

pv = readRDS("~/Desktop/PALVEN_sim_hybclass")
pv = pv[complete.cases(pv$infer2), ]
av = readRDS("~/Desktop/ASPVEN_sim_hybclass")

a = ggplot(pv, aes(x = Var1, y = infer2)) +
  geom_point(aes(size = Freq), col="navy", alpha=.7) +
  scale_size_continuous("Number of individuals",
                        breaks = c(1, 5, 10, 25, 50),
                        range = c(0,10)) +
  xlab("simulated") + ylab("inferred")
b = ggplot(av, aes(x = Var1, y = infer2)) +
  geom_point(aes(size = Freq), col="navy", alpha=.7) +
  scale_size_continuous("Number of individuals",
                        breaks = c(1, 5, 10, 25, 50),
                        range = c(0,10)) +
  xlab("simulated") + ylab("inferred")
ab = (a / b) +  plot_annotation(tag_levels = 'A')
save_plot("~/Desktop/simulations_hybridclass.png", 
          ab, nrow = 2, base_width = 8)


