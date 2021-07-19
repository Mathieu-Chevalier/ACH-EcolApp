formatData <- function(dat, env1=env, poly.range=SP.data, names.sp=NULL, mask=NULL){
  
  ### CH have been used in the two cases but these envelopes are highly sensitive to outliers
  ### We here compute the CH in both spaces as classically done but also a KDE for the geographical space and a MVE for the 
  ### environmental space
  ### We use three measures of distance: euclidean distance to the centroid and to the margins as well as mahalanobis distance 
  
  ### Get information for the species of interest regarding occurence data ###
  #---------------------------------------------------------------------------
  
  # select the species
  abd.data <- filter(dat,Scientific_name==names.sp)
  
  # deal with cases where there are duplicates on latitude and longitude
  abd.data <- abd.data %>% dplyr::group_by(longitude,latitude) %>% dplyr::summarize(scaledAbundance=mean(scaledAbundance,na.rm=TRUE))
  abd.data <- as.data.frame(abd.data)

  # Remove duplicates and thin the data (i.e. points falling within the same cell)
  tmp.dup <- as.data.frame(raster::extract(env1, data.frame(abd.data[,c('longitude','latitude')]), cellnumbers = T))
  ID.cell <- tmp.dup[,1]
  dup <- which(duplicated(ID.cell))

  #---------------
  if(length(dup)!=0){
    tmp.geo.data <- cbind(abd.data, ID.cell)
    tmp.geo.data$ID.cell <- factor(tmp.geo.data$ID.cell, levels=unique(tmp.geo.data$ID.cell))
    tmp.geo.data <- tmp.geo.data %>% dplyr::group_by(ID.cell) %>% 
      dplyr::summarize(scaledAbundance=mean(scaledAbundance, na.rm=TRUE)) # longitude=mean(longitude,na.rm=TRUE), latitude=mean(latitude,na.rm=TRUE),
    tmp.geo.data <- cbind(abd.data[-dup,-3], tmp.geo.data[,-1])
  } else tmp.geo.data <- abd.data
  #---------------
  
  tmp1 <- as.data.frame(tmp.geo.data)
  colnames(tmp1)[1:2] <- c('x','y')
  geo.data <- occ.desaggregation(tmp1, 0.1666667) # this is the resolution of the raster
  colnames(geo.data)[1:2] <- c('longitude','latitude')
  
  # Gather environmental info at sampling plots and remove missing values (e.g. points falling in the ocean)
  env.data <- as.data.frame(raster::extract(env1, data.frame(geo.data[,c('longitude','latitude')])))
  if(any(is.na(env.data))){
    rm1 <- which(is.na(env.data[,1]))
    env.data <- env.data[-rm1,]
    geo.data <- geo.data[-rm1,]
  }
  
  # Only focus on species for which there are at least 10 populations
  if(nrow(na.omit(geo.data)) >= 10){ 
    
    ### Get information over species ranges ###
    #------------------------------------------
    
    # Get the polygon describing the whole species range
    pol <- poly.range[poly.range$binomial==names.sp,]
    line <- as(pol, "SpatialLinesDataFrame")
    grd.pol <- rasterize(pol[,1], mask, 1, background=0, fun='last')
    grd.line <- rasterize(line[,1], mask, 1, background=0, fun='last')
    grd <- max(grd.pol, grd.line)
    
    # Get coordinates over the species range 
    tmp.rg <- as.data.frame(raster::extract(grd, coordinates(grd)))
    tmp.rg <- cbind(coordinates(grd), tmp.rg)
    tmp.rg <- tmp.rg[-which(tmp.rg[,3]==0), 1:2]
    tmp1 <- tmp.rg
    colnames(tmp1)[1:2] <- c('x','y')
    tmp <- occ.desaggregation(tmp1, 0.1666667) # this is the resolution of the raster
    colnames(tmp)[1:2] <- c('longitude','latitude')
    range.data <- tmp
    
    # Extract PCA values over the whole species range
    range.env <- as.data.frame(raster::extract(env1, range.data[,1:2]))
    if(any(is.na(range.env))){
      rm2 <- which(is.na(range.env[,1]))
      range.env <- range.env[-rm2,]
    }
    
    ### Define the different envelopes and compute distances ###
    #-----------------------------------------------------------
    
    ### In the geographic space ###
    
    geo.dist <- NULL
    
    #-------------- CH
    
    # 1- Derive the CH envelope and compute the centroid
    ch <- chull(range.data)
    tmpHull <- range.data[ch,] 
    pss <- SpatialPolygons(list(Polygons(list(Polygon(tmpHull)), ID=1)))
    tmp.pss <- SpatialPolygons2PolySet(pss)
    centroid <- calcCentroid(tmp.pss,rollup=2)[,3:4]
    ID.poly <- rep(1, nrow(geo.data))
    
    # 1.1. Haversine distance relative to centroids
    Dist.ch.cent <- numeric()
    for(z in 1:nrow(geo.data)){Dist.ch.cent[z] <- geosphere::distHaversine(geo.data[z,], centroid)}
    Dist.ch.cent <- Dist.ch.cent/max(Dist.ch.cent, na.rm=T)
    Dist.ch.cent <- cbind(geo.data, env.data, Distance=Dist.ch.cent, space="Geo", envelope="CH", Dist.type="centroid")
    
    # 1.2. Euclidean distance relative to margins
    Dist.ch.marg <- dist.margin(SpatialPoints(geo.data[,1:2]), SpatialPoints(geo.data[,1:2]), niche.delim=pss, ID=ID.poly)$rel.dist
    Dist.ch.marg <- cbind(geo.data, env.data, Distance=Dist.ch.marg, space="Geo", envelope="CH", Dist.type="margin")
    
    # 1.3. Mahalanobis distance relative to centroid
    Sx <- cov(geo.data[,1:2])
    Dist.ch.maha <- mahalanobis(geo.data[,1:2], as.matrix(centroid), Sx)
    Dist.ch.maha <- Dist.ch.maha/max(Dist.ch.maha, na.rm=T)
    Dist.ch.maha <- cbind(geo.data, env.data, Distance=Dist.ch.maha, space="Geo", envelope="CH", Dist.type="mahalanobis")
    
    #-------------- KDE
    
    # 2- Derive the KDE envelope and compute the centroid
    fhat <- kde(x=range.data, compute.cont=TRUE)
    c99 <- contourLines(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, level=fhat$cont[length(fhat$cont)]) 
    l99 <- list()
    ID <- 1:length(c99)
    for(k in 1:length(c99))l99[[k]] <- Polygons(list(Polygon(c99[[k]][-1])), ID=ID[k])
    pss <- gUnaryUnion(SpatialPolygons(l99))
    tmp.pss <- SpatialPolygons2PolySet(pss)
    centroid <- calcCentroid(tmp.pss,rollup=2)[,3:4]
    tmp <- data.frame(ID=1:nrow(range.data))
    points <- SpatialPointsDataFrame(range.data, tmp)
    ID.poly <- raster::extract(pss, points)[,2]
    
    # 2.1. Haversine distance relative to centroids
    Dist.kde.cent <- numeric()
    for(z in 1:nrow(geo.data)){Dist.kde.cent[z] <- geosphere::distHaversine(geo.data[z,], centroid[ID.poly[z],])}
    Dist.kde.cent <- Dist.kde.cent/max(Dist.kde.cent, na.rm=T)
    Dist.kde.cent <- cbind(geo.data, env.data, Distance=Dist.kde.cent, space="Geo", envelope="KDE", Dist.type="centroid")
    
    # 2.2. Euclidean distance relative to margins
    Dist.kde.marg <- dist.margin(SpatialPoints(geo.data[,1:2]), SpatialPoints(geo.data[,1:2]), niche.delim=pss, ID=ID.poly)$rel.dist
    Dist.kde.marg <- cbind(geo.data, env.data, Distance=Dist.kde.marg, space="Geo", envelope="KDE", Dist.type="margin")
    
    # 2.3. Mahalanobis distance relative to centroid
    Sx <- cov(geo.data[,1:2])
    tmp <- data.frame(ID=1:nrow(geo.data))
    pts <- SpatialPointsDataFrame(geo.data[,1:2], tmp)
    tmp.ID <- raster::extract(pss, pts)[,2] 
    
    # Deal with cases where points fall in different polygons
    if(nlevels(factor(tmp.ID))>1){
      Dist.kde.maha <- rep(NA, nrow(geo.data))
      lev <- levels(factor(tmp.ID))
      for(z in 1:length(lev)){
        tmp <- which(tmp.ID==lev[z])
        tmp.dist <- mahalanobis(geo.data[tmp,1:2], as.matrix(centroid[lev[z],]), Sx)
        Dist.kde.maha[tmp] <- tmp.dist
      } 
      if(length(Dist.kde.maha)!=nrow(geo.data)){
        m1 <- which(is.na(Dist.kde.maha))
        m2 <- which(is.na(tmp.ID))
        mat <- match(m1, m2)
        Dist.kde.maha <- Dist.kde.maha[-which(is.na(mat))]
      } 
    } else Dist.kde.maha <- mahalanobis(geo.data[, 1:2], as.matrix(centroid[1,]), Sx)
    Dist.kde.maha <- Dist.kde.maha/max(Dist.kde.maha, na.rm=T)
    Dist.kde.maha <- cbind(geo.data, env.data, Distance=Dist.kde.maha, space="Geo", envelope="KDE", Dist.type="mahalanobis")
    
    # 3. Put all distances in a dataframe
    geo.dist <- rbind(Dist.ch.cent, Dist.ch.marg, Dist.ch.maha, Dist.kde.cent, Dist.kde.marg, Dist.kde.maha)
    
    ### In the environmental space ###
    
    env.dist <- NULL
    
    #-------------- CH
    
    # 1- Derive the CH envelope and compute the centroid
    ch <- chull(range.env)
    tmpHull <- range.env[ch,] 
    pss <- SpatialPolygons(list(Polygons(list(Polygon(tmpHull)), ID=1))) 
    tmp.pss <- SpatialPolygons2PolySet(pss)
    centroid <- calcCentroid(tmp.pss,rollup=2)[,3:4]
    ID.poly <- rep(1, nrow(env.data))
    
    # 1.1. Haversine distance relative to centroids
    Dist.ch.cent <- numeric()
    for(z in 1:nrow(env.data)){Dist.ch.cent[z] <- geosphere::distHaversine(env.data[z,], centroid)}
    Dist.ch.cent <- Dist.ch.cent/max(Dist.ch.cent, na.rm=T)
    Dist.ch.cent <- cbind(geo.data, env.data, Distance=Dist.ch.cent, space="Env", envelope="CH", Dist.type="centroid")
    
    # 1.2. Euclidean distance relative to margins
    Dist.ch.marg <- dist.margin(SpatialPoints(env.data[,1:2]), SpatialPoints(env.data[,1:2]), niche.delim=pss, ID=ID.poly)$rel.dist
    Dist.ch.marg <- cbind(geo.data, env.data, Distance=Dist.ch.marg, space="Env", envelope="CH", Dist.type="margin")
    
    # 1.3. Mahalanobis distance relative to centroid
    Sx <- cov(env.data[,1:2])
    Dist.ch.maha <- mahalanobis(env.data[,1:2], as.matrix(centroid), Sx)
    Dist.ch.maha <- Dist.ch.maha/max(Dist.ch.maha, na.rm=T)
    Dist.ch.maha <- cbind(geo.data, env.data, Distance=Dist.ch.maha, space="Env", envelope="CH", Dist.type="mahalanobis")
    
    #-------------- MVE
    
    # 2- Derive the CH envelope and compute the centroid
    pss <- tolEllipsePlot(range.env)
    tmp.pss <- SpatialPolygons2PolySet(pss)
    centroid <- calcCentroid(tmp.pss,rollup=2)[,3:4]
    ID.poly <- rep(1,nrow(env.data))
    
    # 2.1. Euclidean distance relative to centroids
    Dist.mve.cent <- numeric()
    for(z in 1:nrow(env.data)){Dist.mve.cent[z] <- euc.dist(env.data[z,], centroid)}
    Dist.mve.cent <- Dist.mve.cent/max(Dist.mve.cent, na.rm=T)
    Dist.mve.cent <- cbind(geo.data, env.data, Distance=Dist.mve.cent, space="Env", envelope="MVE", Dist.type="centroid")
    
    # 2.2. Euclidean distance relative to margins
    Dist.mve.marg <- dist.margin(SpatialPoints(env.data[,1:2]), SpatialPoints(env.data[,1:2]), niche.delim=pss, ID=ID.poly)$rel.dist
    Dist.mve.marg <- cbind(geo.data, env.data, Distance=Dist.mve.marg, space="Env", envelope="MVE", Dist.type="margin")
    
    # 2.3. Mahalanobis distance relative to centroid
    Sx <- cov(env.data[,1:2])
    Dist.mve.maha <- mahalanobis(env.data[,1:2], as.matrix(centroid), Sx)
    Dist.mve.maha <- Dist.mve.maha/max(Dist.mve.maha, na.rm=T)
    Dist.mve.maha <- cbind(geo.data, env.data, Distance=Dist.mve.maha, space="Env", envelope="MVE", Dist.type="mahalanobis")
    
    # 3. Put all distances in a dataframe
    env.dist <- rbind(Dist.ch.cent, Dist.ch.marg, Dist.ch.maha, Dist.mve.cent, Dist.mve.marg, Dist.mve.maha)
    
    # Store results
    Out.tab <- data.frame(rbind(geo.dist, env.dist))
    Out.tab$Species <- names.sp
    Out.tab[,1] <- as.numeric(as.character(Out.tab[,1]))
    Out.tab[,2] <- as.numeric(as.character(Out.tab[,2]))
    rownames(Out.tab) <- NULL
  } else Out.tab <- NULL
  
  return(Out.tab)
}



