# This script runs Normal Mixture Models to use unsupervised and supervised
# classification of phenotypic data on 2 hybrid zones of Encelia.
#
# This code is part of the manuscript:
# 
# 	Di Vittorio C, et al. Natural selection mantains species despite widespread 
#		hybridization in the deser shrub Encelia
#
# By: the authors
# Date: June 2020 

# load libraries
library(mclust)
library(tidyverse)

# Load data
indata = read_csv("~/Dropbox/Encelia/manuscripts/EnceliaTransplant_v2/scripts/hybzone_inds_hybclass.csv")

##########
# PALVEN #
##########

# Wrangle data
palven = indata %>%
		filter(hybzone ==  "PALVEN")

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
  dplyr::select(vals) %>% summarize_all(mean, na.rm = T) %>% ungroup()

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
          
# create object for classification analysis
  
tbl = left_join( palven, tbl, 
				 by = c( "ind" = "sample2" ) )
tbl = tbl %>% drop_na()
tbl = tbl %>%
		arrange( species )

# Exploratory Plot
ggplot( tbl, 
	    aes( x = PC1, 
		  color = species) ) +
  	   geom_histogram( fill = "white", 
  			      alpha = 0.5, 
  			      position = "identity")
 
##############
# Unsupervised classification

# Fit NMM for 4 groups on PC1
#mclust.options() # Check Mclust options. For details, see ?mclust.options 
# Save default
opt_mc = mclust.options()
# Use original variables - no transfroamtion 
mclust.options( hcUse = "VARS" )

tbl_nmm = tbl  %>%
	dplyr::select( PC1 ) %>%
	Mclust( G = 1:4 )

summary(tbl_nmm)
tbl_nmm$BIC
plot(tbl_nmm, "BIC")

# Define objects to make other analysis quicker
# Data	
PC1 = pull(tbl, PC1)
# Calssification based on genomic data
Classification = pull(tbl, species)

# Classification table compared to groups assigned by NMMs
table( Classification, 
	   tbl_nmm$classification )
	   
# Adjusted Rand Index: measures of agreement between the two partitions
# a) the genomic classication vs b) the statistical classification 
# based on NMM. 0 = random partition, 1 = perfect accuracy
adjustedRandIndex( Classification, 
				   tbl_nmm$classification )

#plot(tbl_nmm, what = "classification")
#plot(tbl_nmm, what = "uncertainty")

##############
# Supervised classification using the Original classification for training.

# Use NMMs for testing new groups groups
test_classification = MclustDA( PC1, 
							    Classification,
							    G = 1:4 )
							    
# results of classification
summary( test_classification )

# Estimate density (for plotting)
density_nmm = densityMclust( PC1, 
							 G = 1:4 )

# Plot density							 
plot( density_nmm, 
	what = "density")
#plot( density_nmm, 
#	what = "density",
#	type = "persp")
br = seq( min( PC1 ), 
		  max( PC1 ), 
		  length = 21 )
plot( density_nmm, 
	  what = "density",
	  data = PC1,
	  breaks = br )

# Make nicer plot with colored historgrams per "true" classificatioon and densities on top

# First need to recalibrate histogramas per group using Classification
h1 = hist( PC1[Classification == "backcross"], 
			breaks = br, 
			plot = FALSE )
h1$density = h1$density*prop.table(table(Classification))[1]

h2 = hist( PC1[Classification == "F1"], 
			breaks = br,
			plot = FALSE ) 
h2$density = h2$density*prop.table(table(Classification))[2] 

h3 = hist( PC1[Classification == "palmeri"], 
		   breaks = br,
		   plot = FALSE )
h3$density = h3$density*prop.table(table(Classification))[3]

h4 = hist( PC1[Classification == "ventorum"], 
		   breaks = br,
		   plot = FALSE )
h4$density = h4$density*prop.table(table(Classification))[4]

# re estimate values for density
x = seq( min( PC1 ) - diff( range( PC1 )) / 10, 
		 max( PC1 ) + diff( range( PC1 )) / 10,
		 length = 200 )
		 
cdens = predict( density_nmm, 
				  x, 
				  what = "cdens" )
cdens = t(apply( cdens,
				  1,
				  function( d ) d*density_nmm$parameters$pro) )
col = adjustcolor( c("dodgerblue2", "green3", "darkorange", "red3"), 
				   alpha = 0.3 )
				   
# Make plot				   
plot( h1, 
	  xlab = "PC1", 
	  freq = FALSE,
	  main = "",
	  border = FALSE,
	  col = col[2],
	  xlim = range( x ),
	  ylim = range( h1$density, 
	  			    h2$density, 
	  			    cdens ) )
plot( h2, 
	  add = TRUE, 
	  freq = FALSE, 
	  border = FALSE, 
	  col = col[3] )
plot( h3,
	  add = TRUE,
	  freq = FALSE,
	  border = FALSE,
	  col = col[1] )
plot( h4,
	  add = TRUE,
	  freq = FALSE,
	  border = FALSE,
	  col = col[4] )
# Add histograms  
matplot( x, 
		 cdens,
		 type = "l",
		 lwd = 1,
		 add = TRUE,
		 lty = 1:3, 
		 col = 1 )
# decorate		 
box()
legend( "topleft", 
	   legend = c( "E. palmeri",
	   		       "Backcross", 
	   		       "F1",
	   		       "E. ventorum" ),
        col = col,
        pch = 19, 
        cex = 0.65 )


##########
# ASPVEN #
##########

# Wrangle data
aspven = indata %>%
		filter(hybzone ==  "ASPVEN")

inds1 = read_csv("~/Dropbox/Encelia/leaf_images/organized/image_to_sample_map.csv")
inds2 = read_csv("~/Dropbox/Encelia/analysis/hybrid_zones/Encelia_Samples - GENERAL.csv")
# leaf data
p1 = read_csv("~/Dropbox/Encelia/analysis/hybrid_zones/leaf_images/leaf_image - mayra.csv")
p2 = left_join(p1, inds1, by = c("FileName" = "PICTURE_NAME"))

# a few inds are doubled, so merge by sample name not picture name
vals = c("Area", "Mean", "Mode", "Perimeter",
         "Major", "Minor", "Angle", "Circ.", 
         "AR", "Round", "Solidty")
p3 = p2 %>% group_by(SAMPLE_ID, LeafID) %>% 
  dplyr::select(vals) %>% summarize_all(mean, na.rm = T) %>% ungroup()
p3 = p3[complete.cases(p3), ]

# do the pca
pca2 = prcomp(p3 %>% dplyr::select(vals), center = T, scale. = T)
p4 = cbind(p3, pca2$x)
p5 = p4 %>% group_by(SAMPLE_ID) %>% 
  dplyr::select("PC1", "PC2") %>% summarize_all(mean, na.rm = T) %>% ungroup()

# read in data & match up to ind names
tbl = read.table("~/Dropbox/Encelia/analysis/hybrid_zones/NGSadmix/ASPVEN.2.qopt",
                 col.names = c("Ancestry1", "Ancestry2"))
inds = read.csv("~/Dropbox/Encelia/analysis/hybrid_zones/ASPVEN_inds.txt", header=F, stringsAsFactors = F)
tbl$sample = inds$V1
tbl = left_join(tbl, p5, by = c("sample" = "SAMPLE_ID"))

# add in lats and longs
tbl$sample2 = gsub("-","", tbl$sample)
tbl = left_join(tbl, 
                inds2 %>% dplyr::select(PLANT_ID, PUTATIVE_SPECIES, LATITUDE, LONGITUDE, NOTES), 
                by = c("sample2" = "PLANT_ID"))

# remove geographic outliers
tbl = tbl[tbl$LONGITUDE < -114.35, ]

# keep main two parts
hill = tbl[tbl$LONGITUDE < -114.45, ]
tbl = tbl[tbl$LONGITUDE > -114.45, ]
# remove another outlier
tbl = tbl[tbl$LATITUDE < 28.888, ]


# create object for classification analysis

tbl = left_join( aspven, tbl, 
				 by = c( "ind" = "sample2" ) )
tbl = tbl %>% drop_na()
tbl = tbl %>%
		arrange( species )

# Exploratory Plot
ggplot( tbl, 
	    aes( x = PC1, 
		  color = species) ) +
  	   geom_histogram( fill = "white", 
  			      alpha = 0.5, 
  			      position = "identity")
 
##############
# Unsupervised classification

# Fit NMM for 4 groups on PC1
#mclust.options() # Check Mclust options. For details, see ?mclust.options 
# Save default
opt_mc = mclust.options()
# Use original variables - no transfroamtion 
mclust.options( hcUse = "VARS" )

tbl_nmm = tbl  %>%
	dplyr::select( PC1 ) %>%
	Mclust( G = 1:4 )

summary(tbl_nmm)
tbl_nmm$BIC
plot(tbl_nmm, "BIC")

# Define objects to make other analysis quicker
# Data	
PC1 = pull(tbl, PC1)
# Calssification based on genomic data
Classification = pull(tbl, species)

# Classification table compared to groups assigned by NMMs
table( Classification, 
	   tbl_nmm$classification )
	   
# Adjusted Rand Index: measures of agreement between the two partitions
# a) the genomic classication vs b) the statistical classification 
# based on NMM. 0 = random partition, 1 = perfect accuracy
adjustedRandIndex( Classification, 
				   tbl_nmm$classification )

#plot(tbl_nmm, what = "classification")
#plot(tbl_nmm, what = "uncertainty")

##############
# Supervised classification using the Original classification for training.

# Use NMMs for testing new groups groups
test_classification = MclustDA( PC1, 
							    Classification,
							    G = 1:4 )
							    
# results of classification
summary( test_classification )

# Estimate density (for plotting)
density_nmm = densityMclust( PC1, 
							 G = 1:4 )

# Plot density							 
plot( density_nmm, 
	what = "density")
#plot( density_nmm, 
#	what = "density",
#	type = "persp")
br = seq( min( PC1 ), 
		  max( PC1 ), 
		  length = 21 )
plot( density_nmm, 
	  what = "density",
	  data = PC1,
	  breaks = br )

# Make nicer plot with colored historgrams per "true" classificatioon and densities on top

# First need to recalibrate histogramas per group using Classification
h1 = hist( PC1[Classification == "asperifolia"], 
			breaks = br, 
			plot = FALSE )
h1$density = h1$density*prop.table(table(Classification))[1]

h2 = hist( PC1[Classification == "backcross"], 
			breaks = br,
			plot = FALSE ) 
h2$density = h2$density*prop.table(table(Classification))[2] 

h3 = hist( PC1[Classification == "F1"], 
		   breaks = br,
		   plot = FALSE )
h3$density = h3$density*prop.table(table(Classification))[3]

h4 = hist( PC1[Classification == "ventorum"], 
		   breaks = br,
		   plot = FALSE )
h4$density = h4$density*prop.table(table(Classification))[4]

# re estimate values for density
x = seq( min( PC1 ) - diff( range( PC1 )) / 10, 
		 max( PC1 ) + diff( range( PC1 )) / 10,
		 length = 200 )
		 
cdens = predict( density_nmm, 
				  x, 
				  what = "cdens" )
cdens = t(apply( cdens,
				  1,
				  function( d ) d*density_nmm$parameters$pro) )
col = adjustcolor( c("dodgerblue2", "green3", "darkorange", "red3"), 
				   alpha = 0.3 )
				   
# Make plot				   
plot( h1, 
	  xlab = "PC1", 
	  freq = FALSE,
	  main = "",
	  border = FALSE,
	  col = col[1],
	  xlim = range( x ),
	  ylim = range( h1$density, 
	  			    h1$density, 
	  			    cdens ) )
plot( h2, 
	  add = TRUE, 
	  freq = FALSE, 
	  border = FALSE, 
	  col = col[2] )
plot( h3,
	  add = TRUE,
	  freq = FALSE,
	  border = FALSE,
	  col = col[3] )
plot( h4,
	  add = TRUE,
	  freq = FALSE,
	  border = FALSE,
	  col = col[4] )
# Add histograms  
matplot( x, 
		 cdens,
		 type = "l",
		 lwd = 1,
		 add = TRUE,
		 lty = 1:3, 
		 col = 1 )
# decorate		 
box()
legend( "topleft", 
	   legend = c( "E. asperifolia",
	   		       "Backcross", 
	   		        "F1",
	   		        "E. ventorum" ),
        col = col,
        pch = 19, 
        cex = 0.65 )
