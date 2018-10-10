#' Create a density  plot of posterior distribution for a given parameter
#' 
#' @description  creates a density plot for a parameter
#' 
#' @param col column that contains the parameter posterior draws
#' @param Mlist list that contains posterior draws from several gibbs sampler runs
#' 
#' @importFrom stats density
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot lines
#' @importFrom stringr str_sub str_replace_all
#' 
#' @export

density_plot=function(col,Mlist){

  dir.create("results/densityplots")
  colname=str_replace_all(col,"[.]","_")
  imagename=paste0("results/densityplots/",colname,".png")
  png(imagename)
  plot(density(Mlist[[1]][,col]),type="l",
       main=loc_hash[[str_sub(col,-2)]],xlab="x",ylab="")
  lines(density(Mlist[[2]][-(1:Nburn),col]),col="red")
  lines(density(Mlist[[3]][-(1:Nburn),col]),col="green")
  dev.off()
}