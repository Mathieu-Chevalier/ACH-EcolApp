extract.slope <- function(x, y, data){
  mod <- lm(y ~ x, data=data)
  return(coef(mod)[2])
}


