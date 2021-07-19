gamR2 <- function(gam){
  R2 <- 1-((sum(residuals(gam)^2))/(sum((gam$y - mean(gam$y))^2)))
  R2adj <- 1- ((1 - R2) * (length(gam$y) - 1)/(length(gam$y) - length(gam$coefficients) - 1))
  a <- data.frame(R2, R2adj)
  return(a)
}