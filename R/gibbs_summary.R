#' Create summary output for Gibbs Sampler
#' 
#' @description  creates a trace plot for posterior draws
#' 
#' @param x column of posterior draws for a parameter
#' @export
#' 
#'
 
gibbs_summary=function(x){ 
  summ=c(min(x),quantile(x,c(0.025,0.25)),median(x),mean(x),sd(x),
         quantile(x,1-rev(c(0.025,0.25))),max(x))
  names(summ)=c("min","2.5%","25%","Median","Mean","sd","75%","97.5%","max")
  summ
}