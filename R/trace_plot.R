#' Create a trace plot of posterior draws for a given parameter
#' 
#' @description  creates a trace plot for posterior draws
#' 
#' @param col column that contains the parameter posterior draws
#' @param Mlist list that contains posterior draws from several gibbs sampler runs
#' 
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot lines
#' @importFrom stringr str_sub str_replace_all
#' 
#' @export

trace_plot=function(col,Mlist){
  
  dir.create("results/traceplots")
  colname=str_replace_all(col,"[.]","_")
  imagename=paste0("results/traceplots/",colname,".png")
  png(imagename)
  plot(Mlist[[1]][,col],type="l",xlab="Iterations",ylab=col)
  lines(Mlist[[2]][,col],col="red")
  lines(Mlist[[3]][,col],col="green")
  dev.off()
}