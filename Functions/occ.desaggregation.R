occ.desaggregation=function (xy, min.dist, by = NULL) {
    if (is.null(xy$x) | is.null(xy$y)) {
        stop("no x and/or y column")
    }
    if (!is.null(by)) {
        xy[, which(names(xy) == by)] <- factor(xy[, which(names(xy) == 
            by)])
    }
    new.data <- NULL
    if (is.null(by)) {
        to <- 1
    }
    else {
        to <- nlevels(xy[, which(names(xy) == by)])
    }
    for (i in 1:to) {
        if (to > 1) {
            print(paste("desaggregate species", i))
        }
        if (!is.null(by)) {
            del.min.dist <- xy[xy[, by] == levels(xy[, by])[i], 
                ]
        }
        else {
            del.min.dist <- xy
        }
        if (sum(duplicated(paste(del.min.dist$x, del.min.dist$y))) > 
            0) {
            stop(paste("duplicated values", levels(xy[, 
                by])[i], sep = " "))
        }
        repeat {
            nn1 <- round(nndist(del.min.dist[, "x"], del.min.dist[, 
                "y"]), nchar(strsplit(as.character(min.dist), "\\.")[[1]][2]))
            if (sum(nn1 < min.dist) == 0) {
                break
            }
            nn2 <- round(nndist(del.min.dist[, "x"], del.min.dist[, 
                "y"], k = 2), nchar(strsplit(as.character(min.dist), "\\.")[[1]][2]))
            del1 <- nn1 == min(nn1)
            del2 <- nn2 == min(nn2[del1])
            delk <- del1 & del2
            if (sum(del2) > 1) {
                for (k in 3:8) {
                  nn <- nndist(del.min.dist[, "x"], del.min.dist[, 
                    "y"], k = k)
                  delk <- delk & nn == min(nn[delk])
                  if (sum(nn[delk] == min(nn[delk])) > 1) {
                    break
                  }
                }
            }
            del.min.dist <- del.min.dist[-(which(delk)[1]), ]
        }
        new.data <- rbind(new.data, del.min.dist)
    }
    result <- list(initial = nrow(xy), kept = nrow(new.data), 
        out = nrow(xy) - nrow(new.data))
    return(xy = new.data[, !colnames(new.data) %in% c("nn", 
        "nn2", "id")])
}