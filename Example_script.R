##############################################################################################################
# This script aims to measure the distance of populations using different methodological settings and to
# subsequently evaluate the distance-abundance relationship under these settings using different methods.
#
# The settings to compute distances involve:
#           - Distance measured relative to centroids (Euclidean and Mahalanobis) or margins
#           - Distance measured within an envelope characterized with a convex hull (CH), a kernel density estimator (KDE) or a minimum volume ellipsoid (MVE)
#           - Distance measured in the geographical or in the environmental space
#
# The methods to characterize the distance abundance relationship involve:
#           - Correlation coefficients measured with Spearman or Pearson
#           - Generalized additive models with random effects on intercepts and splines
#           - Bayesian hierarchical model with random effects on intercepts, slopes and quadratic terms
#
# This script is applied to the mammal dataset and can be run effectively within a reasonable amount of time.
# However, for larger datasets (e.g. both bird datasets and the tree dataset), computation time can be very prohibitive and we recommend 
# to use parallel computing in this case.
##############################################################################################################

rm(list=ls());gc() # clean working environment

### load libraries
require(raster)
require(ade4)
require(data.table)
require(dplyr)
require(boa)
require(rgeos)
require(ks)
require(tidymv)
require(ggplot2)
require(robustbase)
require(car)
require(spatstat)
require(ecospat)
require(reshape2)
require(mgcv)
require(gratia)
require(runjags)
require(R2jags)
require(qpcR)
require(plyr)
require(ggpubr)
require(raptr)
require(PBSmapping)

### Define workdir
PATH <- "Z:/projects-unil/DALLAS/Files_to_submit/"
setwd(PATH)

### Load some home-made functions
source("Functions/dist.margin.R")
source("Functions/formatData.R")
source("Functions/tolEllipsePlot.R")
source("Functions/occ.desaggregation.R")
source("Functions/euc.dist.R")
source("Functions/gamR2.R")

### Define the group for which inferences are to be drawn
all.group <- c("mammals", "birds", "fishes", "trees", "BBSbirds")
focus.group <- all.group[1] #### Mammals

### Create a folder to store outputs
if(!file.exists(paste("outputs", focus.group, sep="/"))) dir.create(paste("outputs", focus.group, sep="/"), recursive = T)

### Define methodological settings
tmp1 <- expand.grid(c("centroid", "margin", "mahalanobis"), c("CH", "KDE"), "Geo")
tmp2 <- expand.grid(c("centroid", "margin", "mahalanobis"), c("CH", "MVE"), "Env")
settings <- rbind(tmp1, tmp2)
name.settings <- apply(settings, 1, function(x)paste(sapply(x, as.character), collapse="_"))

### Read abundance data for the considered group
data <- readRDS(paste0("Data/Abundances/data_", focus.group, ".Rdata"))

### Species names with at least 10 occurences
tmp <- table(data$Scientific_name)
names.10.occ <- names(which(tmp >= 10))

#####################################################
### 0.1 Get PCA information and mask
#####################################################

mask <- raster("Data/mask_10min.tif") # mask with a 10min resolution
env.stack <- readRDS("Data/PCA.env.Rdata") # PCA data resulting from an analysis performed on 19 climatic layers extracted from worldclim at a 10 min resolution 

#####################################################
### 0.2 Load IUCN range maps for the different species
#####################################################

### Range maps for trees were downloaded from 
### https://www.fs.fed.us/nrs/atlas/littlefia/species_table.html
### https://www.fs.fed.us/nrs/atlas/littlefia/

SP.data <- get(load(paste0("Data/Shapefiles/", focus.group, ".Rdata")))
if(any(colnames(SP.data@data)=="SCINAME")) colnames(SP.data@data)[which(colnames(SP.data@data)=="SCINAME")] <- "binomial"
names.range <- levels(factor(SP.data$binomial))

### Get common names between abundance data and range maps
common.names <- intersect(names.10.occ, names.range)
data <- data[which(data$Scientific_name %in% common.names), ] # Only keep species in common between the two datasets

############################################
### 1. Compute distances for the given group
############################################

file.create(paste0("outputs/", focus.group, "/out.distances.txt")) # create empty file to store results
for(i in 1:length(common.names)){
  tmp.out <- formatData(dat=data, env1=env.stack, poly.range=SP.data, names.sp=common.names[i], mask=mask) 
  if(!is.null(tmp.out)) write.table(tmp.out, file=paste0("outputs/", focus.group, "/out.distances.txt"), append=T, col.names=FALSE)
  print(length(common.names)-i)
}

############################################
### 2. Clean distance data (e.g. remove species not having 10 occurences) 
############################################

df <- read.table(paste0("outputs/", focus.group, "/out.distances.txt"), row.names=NULL, fill=T)
df <- df[,-1]
colnames(df) <- c("longitude", "latitude", "scaledAbundance", "Axis1", "Axis2", "Distance", "space", "envelope", "Dist.type", "Species")

### Remove some values that are no longer of interest (only useful to compute species envelopes)
if(any(na.omit(df$scaledAbundance)==0)) df <- df[-which(df$scaledAbundance==0),] # Remove zero values in scaledAbundances 
if(any(is.na(df$scaledAbundance))) df <- df[-which(is.na(df$scaledAbundance)),] # Remove missing values in scaledAbundances

### remove populations for which we haven't been able to compute distance for some methodological settings -> homogeneize populations across datasets
if(any(is.na(df$Distance))){
  a  <- df[which(is.na(df$Distance)),]
  rm <- unique(paste(a[,1], a[,2], a[,10], sep="_"))
  df$pasting <- paste(df[,1], df[,2], df[,10],sep="_")
  df <- df[-which(df$pasting %in% rm),]
  df <- df[,-ncol(df)]
}

### Remove species not having at least 10 occurrences
sp.occ <- data.frame(df %>% dplyr::group_by(space, envelope, Dist.type, Species) %>% dplyr::summarize(Distance = n()))
tmp <- unique(sp.occ[,4:5])
names.rm <- levels(factor(tmp[which(tmp[,2]<10),1]))
df <- df[-which(df$Species %in% names.rm),]
df$Species <- factor(df$Species)

### Save cleaned data
saveRDS(df, file=paste0("outputs/", focus.group, "/out.distances_clean.Rdata"))

################################################################################
### 3. Compute correlations between distance and abundance (Spearman)
################################################################################

df <- readRDS(paste0("outputs/", focus.group, "/out.distances_clean.Rdata"))

### Spearman correlations
Cor.spearman <- df %>% dplyr::group_by(space, envelope, Dist.type, Species) %>% dplyr::summarize(Correl=cor(scaledAbundance, Distance, method="spearman"))

### Evaluate whether correlation coefficients are different from zero
stat.test <- compare_means(Correl~1, mu=0, Cor.spearman, group.by=c("space", "envelope", "Dist.type"))

### The plot
plot.correl <- ggplot(Cor.spearman, aes(x=space, y=Correl, color=envelope))+ 
  geom_violin(position=position_dodge(width=1), aes(fill=envelope), alpha=0.2, trim=FALSE)+
  stat_summary(fun.data=mean_sdl, geom="pointrange", position=position_dodge(width=1))+
  theme_bw()+
  labs(shape="Space", colour="Envelope", y="Spearman's correlations", x="")+
  facet_wrap(~Dist.type)+
  geom_text(data = stat.test, aes(y=1,label = p.signif), position=position_dodge(width=1))+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_manual(values=c("#de2d26", "#31a354", "#2b8cbe"), name="Enveloppe", labels = c("CH", "KDE", "MVE")) + 
  scale_fill_manual(values=c("#de2d26", "#31a354", "#2b8cbe"), name="Enveloppe", labels= c("CH", "KDE", "MVE"))+
  coord_fixed(ratio=0.8)

##########################################################################
### 4 Investigate the distance abundance relationship using LMER with GAMs
##########################################################################

### Get the data
df <- readRDS(paste0("outputs/", focus.group, "/out.distances_clean.Rdata"))
df$scaledAbundance <- log(df$scaledAbundance) # Log abundances

#####################################
# Draw mixed effect models using GAMs
#####################################

for(i in 1:nrow(settings)){
  tmp <- df[which(df$Dist.type==settings[i,1] & df$envelope==settings[i,2] & df$space==settings[i,3]),]
  mod.gam <- bam(scaledAbundance ~ s(Distance, k=3) + s(Distance, Species, bs="fs", k=3), data=tmp, method="REML")
  saveRDS(mod.gam, file=paste0("outputs/", focus.group, "/GAM_",name.settings[i],".Rdata"))
  print(nrow(settings) - i)
}

#####################################
### Compute evaluation metrics
#####################################

mod.gam <- list()
for(i in 1:nrow(settings)) mod.gam[[i]] <- readRDS(paste0("outputs/", focus.group, "/GAM_",name.settings[i],".Rdata"))
names(mod.gam) <- name.settings

### R²
tmp.Rsq <- round(unlist(lapply(mod.gam, function(x)gamR2(x)[1])), 2)
tmp.Rsq.adj <- round(unlist(lapply(mod.gam, function(x)gamR2(x)[2])), 2)

### AIC, AIC weights and delta AIC
tmp.aic <- unlist(lapply(mod.gam, function(x)AIC(x)))
pos <- order(tmp.aic)
tmp.aic.wt <- akaike.weights(tmp.aic)
tmp.aic.df <- cbind(tmp.aic, do.call("cbind", tmp.aic.wt))[,-3]
tmp.aic.df <- apply(tmp.aic.df, 2, function(x) formatC(x, format="f", digits=2))
df.out.gam <- cbind(name.settings, tmp.Rsq, tmp.Rsq.adj, tmp.aic.df)
df.out.gam <- df.out.gam[pos,]
rownames(df.out.gam) <- NULL
colnames(df.out.gam) <- c("Model name", "R²", "R².adj", "AICc", "delta AICc", "AICc weights")
df.out.gam <- as.data.frame(df.out.gam)

### Save results
write.table(df.out.gam, paste0("outputs/", focus.group, "/Summary_GAM_models.txt"))

#####################################
### Perform species-wise predictions 
#####################################

### The predictions
Preds <- lapply(mod.gam, function(x)as.data.frame(predict_gam(x)))
df.sp <- rbindlist(Preds)
vec.names <- strsplit(sapply(names(Preds), function(x)rep(x, nrow(Preds[[1]]))), "_")
df.sp$Dist.type <- unlist(lapply(vec.names, function(x)x[1]))
df.sp$envelope <- unlist(lapply(vec.names, function(x)x[2]))
df.sp$space <- unlist(lapply(vec.names, function(x)x[3]))
df.sp$methods <- paste(df.sp$envelope, df.sp$space, df.sp$Dist.type)
saveRDS(df.sp, paste0("outputs/", focus.group, "/PredGAM_sp_wise.Rdata"))

### Extract slope coefficient from a linear model to add an indicative color scale to the figure to highlight the strength of the relationship
df.pred <- readRDS(paste0("outputs/", focus.group, "/PredGAM_sp_wise.Rdata"))
df.pred$temp <- as.factor(paste(df.pred$methods, df.pred$Species))
tmp.coef <- with(df.pred,by(df.pred, temp, function(x) lm(fit ~ Distance, data = x)))
coef <- sapply(tmp.coef, function(x)coef(x)[2])
names(coef) <- gsub(".Distance", "", names(coef))
DF <- data.frame(cbind(coef=as.numeric(coef), temp=names(coef)))
DF$coef <- as.numeric(DF$coef)

### The plot
df.pred2 <- left_join(df.pred, DF)
sp.plot <- ggplot(data=df.pred2, aes(x=Distance, y=fit, color=coef, group=Species, alpha=abs(coef)))+
  geom_line(size=0.5)+
  facet_wrap(~methods, ncol=3)+
  theme_bw()+
  labs(y="Log(Predicted abundance)", x="Standardized distance", color="Coefficients")+
  scale_alpha_continuous(range=c(0.6,1))+
  scale_colour_gradient2(mid="grey60", low="red", high="green")+
  guides(alpha = "none")+
  coord_fixed(ratio=0.15)

#####################################
### Perform global predictions 
#####################################

### The predictions
Preds.1 <- lapply(mod.gam, function(x)as.data.frame(predict_gam(x, exclude_terms = "s(Distance, Species)", values=list(Species=NULL))))
Preds.2 <- lapply(Preds.1, function(x){
  x$conf.low <- x$fit-1.96*x$se
  x$conf.high <- x$fit+1.96*x$se
  return(x)
})
df.glob <- rbindlist(Preds.2)
vec.names <- strsplit(sapply(names(Preds.2), function(x)rep(x, nrow(Preds.2[[1]]))), "_")
df.glob$Dist.type <- unlist(lapply(vec.names, function(x)x[1]))
df.glob$envelope <- unlist(lapply(vec.names, function(x)x[2]))
df.glob$space <- unlist(lapply(vec.names,function(x)x[3]))
saveRDS(df.glob, paste0("outputs/", focus.group, "/PredGAM_global.Rdata"))

### The plot
df.glob <- readRDS(paste0("outputs/", focus.group, "/PredGAM_global.Rdata"))
df.glob$method <- paste(df.glob$space, df.glob$Dist.type)

global.plot <- ggplot(df.glob, aes(x=Distance, y=fit, color=envelope, fill=envelope))+
  geom_line(size=1.5)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.3, linetype=0)+
  scale_colour_manual(values=c("#de2d26", "#31a354", "#2b8cbe"))+
  scale_fill_manual(values=c("#de2d26", "#31a354","#2b8cbe"))+
  labs(fill="Envelope", colour="Envelope")+
  theme_bw()+
  facet_wrap(~method, ncol=6)+
  labs(y="Log(Predicted abundance)", x="Standardized distance")+
  coord_fixed(ratio=0.4)

#################################################################################
### 5. Investigate the distance abundance relationship using a multi-level 
### modelling framework with random effects on both intercepts and slopes (JAGS)
#################################################################################

### Get the data
df <- readRDS(paste0("outputs/", focus.group, "/out.distances_clean.Rdata"))
df$scaledAbundance <- log(df$scaledAbundance) # Log transform abundances

### Define MCMC settings
store <- 1000
nadap <- 1000
chains <- 3
nburn <- 5000
thin <- 10 

### Variables to monitor
my_var <- c("alpha", "beta", "mean.alpha", "mean.beta", "rho", "mean.quad", "beta.quad", "sd.quad", "sd.alpha", "sd.beta", "sigma", "bpvalue")

#####################################
### Run the models
#####################################

sp <- levels(factor(df$Species))
Nsp <- length(sp)

for(i in 1:nrow(settings)){
  
  tmp <- df[which(df$Dist.type==settings[i,1] & df$envelope==settings[i,2] & df$space==settings[i,3]),]
  
  ### prepare the data to feed JAGS 
  Npop <- table(tmp$Species)
  Dist <- Y <- array(dim=c(Nsp, max(Npop)))
  for(l in 1:Nsp){
    A <- subset(tmp, Species==sp[l])
    Dist[l, 1:Npop[l]] <- A$Distance
    Y[l, 1:Npop[l]] <- A$scaledAbundance
  }
  data <- list(Dist=Dist, Npop=as.vector(Npop), Nsp=Nsp, Y=Y)
  
  ### Run the model
  out0 <- jags.parallel(data=data, parameters.to.save=my_var, model.file="JAGS_models/LMER_quadratic.bug",
                        n.chains=chains, n.burnin=nburn, n.iter=(store*thin)+nburn, n.thin=thin, 
                        export_obj_names = c("data", "nburn", "store", "thin"))
  
  ### Save model outputs
  saveRDS(out0, file=paste0("outputs/", focus.group, "/JAGS_", name.settings[i], ".Rdata"))
  
  print(nrow(settings) - i)
  
}

##############################
### Compute evaluation metrics
##############################

df.out <- NULL
df <- readRDS(paste0("outputs/", focus.group, "/out.distances_clean.Rdata"))

for(i in 1:nrow(settings)){
  
  out0 <- readRDS(paste0("outputs/", focus.group, "/JAGS_", name.settings[i], ".Rdata"))
  out1 <- as.mcmc(out0)
  comb <- combine.mcmc(out1)
  
  ### Get Bayesian p-value
  Bpval <- mean(comb[,"bpvalue"])
  
  ### Get estimated coefficients and associated 95% HPD intervals
  median.slp <- round(median(comb[,"mean.beta"]), 2) 
  HDI.slp <- as.vector(round(boa.hpd(comb[,"mean.beta"], alpha=0.05), 2))
  pst.slp <- paste(median.slp, "[",HDI.slp[1],";", HDI.slp[2], "]" ,sep="")
  median.quad <- round(median(comb[,"mean.quad"]), 2)
  HDI.quad <- as.vector(round(boa.hpd(comb[,"mean.quad"], alpha=0.05), 2))
  pst.quad <- paste(median.quad, "[", HDI.quad[1], ";", HDI.quad[2], "]", sep="")
  
  ### Compute posterior probabilities of coefficients being different from zero
  beta <- comb[,"mean.beta"]
  quad <- comb[,"mean.quad"]
  if(median.slp>0) perc.post.slp <- length(which(beta>0))/length(beta)*100 else perc.post.slp <- length(which(beta<0))/length(beta)*100
  if(median.quad>0) perc.post.quad <- length(which(quad>0))/length(quad)*100 else perc.post.quad <- length(which(quad<0))/length(quad)*100
  
  ### Compute posterior probabilities for expected effects 
  if(median.slp>0){
    perc.neg <- length(which(quad < -(beta)))/length(beta)*100
  } else {
    perc.neg <- length(which(quad < abs(beta)))/length(beta)*100
  } 
  perc.pos <- 100-perc.neg
  if(settings[i, 1] == "margin") perc <- perc.pos else perc <- perc.neg 
  
  ### Compute marginal and conditional R-square
  tmp <- df[which(df$Dist.type==settings[i,1] & df$envelope==settings[i,2] & df$space==settings[i,3]),]
  yPredFixed <- median(comb[,"mean.alpha"]) + median.slp*tmp$Distance + median.quad*(tmp$Distance^2)
  varFixed <- var(yPredFixed) # variance explained by fixed effects
  varResidual <- median(comb[,"sigma"]^2) # Residual variance
  varRandom <- median(comb[,"sd.alpha"]^2) + median(comb[,"sd.beta"]^2)  # variance explained by random effects
  
  marginalR2 <- varFixed / (varFixed + varRandom + varResidual)
  conditionalR2 <- (varRandom + varFixed) / (varFixed + varRandom + varResidual) 
  
  med.marginal <- round(median(marginalR2), 4)
  med.conditional <- round(median(conditionalR2), 4)
  
  ### Store results
  df.out <- rbind(df.out, c(Dist.type=as.character(settings[i,1]), envelope=as.character(settings[i,2]), 
                            space=as.character(settings[i,3]), group=all.group[1],
                            Bpval, med.marginal, med.conditional, pst.slp, perc.post.slp, pst.quad, perc.post.quad, perc))
  
  print(nrow(settings)-i)
  
}

df.out <- as.data.frame(df.out)
colnames(df.out) <- c("Dist.type", "Space", "Envelope", "Group", "Bayesian p-value", "marginalR2", "conditionalR2",
                      "Distance effect","Post.prob.slope","Quadratic effect","Post.prob.quad", "Pexpected.effect")
write.table(df.out, paste0("outputs/", focus.group, "/Summary_JAGS_models.txt"))

##############################
### Perform species-wise predictions 
##############################

### Get data
df <- readRDS(paste0("outputs/", focus.group, "/out.distances_clean.Rdata"))
Nsp <- nlevels(factor(df$Species))
sp <- levels(factor(df$Species))
file.create(paste0("outputs/", focus.group, "/PredJAGS_sp_wise.txt")) # create empty file to store results

for(i in 1:nrow(settings)){
  
  ### Extract distance data
  tmp <- df[which(df$Dist.type==settings[i,1] & df$envelope==settings[i,2] & df$space==settings[i,3]),]
  SEQ <- seq(min(tmp$Distance), max(tmp$Distance), length.out=100)
  
  ### Extract posterior distributions
  out0 <- readRDS(paste0("outputs/", focus.group, "/JAGS_", name.settings[i],".Rdata"))
  out0 <- as.mcmc(out0) # transform outputs to MCMC objects
  comb <- combine.mcmc(out0) # combine MCMC chains
  
  ### Compute predictions and support for expected relationship
  for(j in 1:Nsp){
    
    if(any(c(i,j)!=1)) name.col <- FALSE else name.col <- TRUE
    
    #-- Predictions
    inter.sp <- median(comb[,paste0("alpha[",j,"]")] + comb[,"mean.alpha"])
    slp.sp <- median(comb[,paste0("beta[",j,"]")] + comb[,"mean.beta"])
    quad.sp <- median(comb[,paste0("beta.quad[",j,"]")] + comb[,"mean.quad"])
    Y.pred <- inter.sp + slp.sp*SEQ + quad.sp*SEQ^2
    
    #-- Support for the ACH
    beta.sp <- comb[,paste0("beta[",j,"]")] + comb[,"mean.beta"]
    quad.sp <- comb[,paste0("beta.quad[",j,"]")] + comb[,"mean.quad"]
    if(slp.sp>0){
      perc.neg <- length(which(quad.sp < -(beta.sp)))/length(beta.sp)*100
    } else {
      perc.neg <- length(which(quad.sp < abs(beta.sp)))/length(beta.sp)*100
    } 
    perc.pos <- 100-perc.neg
    if(settings[i, 1] == "margin") perc <- perc.pos else perc <- perc.neg 
    
    #-- Store results
    Predicted <- cbind(Y=Y.pred, dist=SEQ, P.exp=perc, species=sp[j], Dist.type=as.character(settings[i,1]),
                       envelope=as.character(settings[i,2]), space=as.character(settings[i,3]))
    write.table(Predicted, file=paste0("outputs/", focus.group, "/PredJAGS_sp_wise.txt"), append=T, col.names=name.col)
  }
  
  print(nrow(settings)-i)
  
}

### The plot
df.pred <- read.table(paste0("outputs/", focus.group, "/PredJAGS_sp_wise.txt"), h=T, row.names = NULL)
df.pred$methods <- paste(df.pred$space, df.pred$envelope, df.pred$Dist.type, sep=" ")

plot.sp <- ggplot(data=df.pred, aes(x=dist, y=Y, color=P.exp, group=species, alpha=P.exp))+
  geom_line(size=0.5)+
  facet_wrap(~methods, ncol=3)+
  theme_bw()+
  labs(x="Standardized distance", y="Log(Predicted abundance)", color="P[expected] (%)")+
  scale_alpha_continuous(range=c(0.2, 1))+
  scale_colour_gradient2(mid="grey60",low="red",high="green", midpoint=50)+
  coord_fixed(ratio=0.15)+
  guides(alpha="none")

##############################
### Perform global predictions 
##############################
 
df <- readRDS(paste0("outputs/", focus.group, "/out.distances_clean.Rdata"))
file.create(paste0("outputs/", focus.group, "/PredJAGS_global.txt")) # create empty file to store results

for(i in 1:nrow(settings)){
  
  ### Extract distance data
  tmp <- df[which(df$Dist.type==settings[i,1] & df$envelope==settings[i,2] & df$space==settings[i,3]),]
  SEQ <- seq(min(tmp$Distance), max(tmp$Distance), length.out=100)
  
  if(i!=1) name.col=FALSE else name.col=TRUE
  
  ### Extract posterior estimates of fixed effects
  out0 <- readRDS(paste0("outputs/", focus.group, "/JAGS_", name.settings[i], ".Rdata"))
  out0 <- as.mcmc(out0)
  comb <- combine.mcmc(out0)
  
  ### Draw the predictions 
  glob.inter <- comb[,"mean.alpha"]
  glob.slp <- comb[,"mean.beta"]
  glob.quad <- comb[,"mean.quad"]
  pred <- glob.inter + outer(glob.slp, SEQ, "*") + outer(glob.quad, SEQ^2, "*") # posterior predictions
  pred.HPD <- apply(pred, 2, function(x) quantile(x, probs=c(0.025, 0.975))) # 95% HPD interval of predictions
  pred.med <- apply(pred, 2, median) # median of posterior predictions
  pred.low <- pred.HPD[1,]
  pred.high <- pred.HPD[2,]
  
  ### Store results
  tmp.pred <- cbind(Y=pred.med, Y.low=pred.low, Y.high=pred.high, Distance=SEQ,
                    Dist.type=as.character(settings[i,1]), envelope=as.character(settings[i,2]), 
                    space=as.character(settings[i,3]), group=focus.group)
  write.table(tmp.pred, file=paste0("outputs/", focus.group, "/PredJAGS_global.txt"), append=T, col.names=name.col)
  
  print(nrow(settings)-i)
}

### The plot
df.pred <- read.table(paste0("outputs/", focus.group, "/PredJAGS_global.txt"), h=T, row.names=NULL)
df.pred$methods <- paste(df.pred$space, df.pred$Dist.type)

plot.glob <- ggplot(df.pred, aes(x=Distance, y=Y, ymin=Y.low, ymax=Y.high, color=envelope, fill=envelope)) +
  geom_line(size=1.5)+
  geom_ribbon(aes(ymin=Y.low, ymax=Y.high), alpha=0.3, linetype=0)+
  scale_colour_manual(values=c("#de2d26", "#31a354", "#2b8cbe"))+
  scale_fill_manual(values=c("#de2d26", "#31a354","#2b8cbe"))+
  theme_bw()+
  facet_wrap(~methods, ncol=6)+
  coord_fixed(ratio=0.4)+
  labs(y="Log(Predicted abundance", x="Standardized distance", color="Envelope", fill="Envelope")

#############################################################################################################################
#### END SCRIPT
#############################################################################################################################
