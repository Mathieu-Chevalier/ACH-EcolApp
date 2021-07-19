dist.margin <- function(foc.pop,sp.dist,niche.delim,ID){

  sp.line <- as(niche.delim, 'SpatialLines')
  foc.pop.inout <- over(foc.pop,niche.delim)
  sp.dist.inout <- over(sp.dist,niche.delim)
  sp.inn <- gDistance(spgeom1=sp.dist, spgeom2=sp.line, byid=T)
  sp.inn[which(is.na(sp.dist.inout))] <- NA 
  foc.inn <- gDistance(spgeom1=foc.pop, spgeom2=sp.line, byid=T)
  foc.inn[which(is.na(foc.pop.inout))] <-  NA 
  
  tmp.foc <- tmp.sp <- numeric()
  for(i in 1:ncol(foc.inn)){
		tmp.foc[i] <- foc.inn[ID[i],i]
		tmp.sp[i] <- sp.inn[ID[i],i]
  }
  foc.inn.rel <- tmp.foc/max(tmp.sp,na.rm=T)

  return(list(rel.dist=foc.inn.rel, abs.dist=tmp.foc))
}

