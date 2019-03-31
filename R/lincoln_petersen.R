#' Calculates the Lincoln-Petersen estimate, variance, and 95% CI
#' 
#' @description Calculates the Lincoln-Petersen estimate, variance, and 95% CI
#' 
#' @param n1 the key pop count from the first survey (or listing) 
#' @param n2 the key pop count from the second survey (or listing)
#' @param m the count of individuals in both the first and second surveys or listings
#' 
#' @return a list with element 1 the lincoln-petersen estimate, element 2 the variance, and element 3 the 95% confidence interval
#' 
#' 
#' @export

lincoln_petersen = function(n1, n2, m) {
  
  n1 = as.numeric(n1)
  n2 = as.numeric(n2)
  m = as.numeric(m)
  #list to contain the estimate, variance, and CI
  return.list = list()
  #Lincoln-Petersen estimator
  return.list[[1]] = (n1 * n2)/m
  #variance of lincoln-petersen estimator
  return.list[[2]] = (n1 * n2 * (n1 - m) * (n2 -m))/m^3
  #normal approximation of 95% CI
  return.list[[3]] = return.list[[1]] + c(-1.96, 1.96) * sqrt(return.list[[2]])
  
  return(return.list)
}