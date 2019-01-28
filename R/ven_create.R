#' Creates a Venn Diagram with either  2 or 3 overlaps
#' 
#' @description Creates a listing Venn Diagram showing 
#' @description missing and known data
#' 
#' @param nListings number of listings
#' @param vecCounts a vector that contains the known number of participants for each partition of the Venn diagram
#' @param listNames a vector of the listing names for the legend
#' @param locName the location name for this Venn Diagram
#' @param groupColors color of the circle for each listing
#' 
#' @return a ggplot object
#' 
#' @import ggplot2
#' @import ggforce
#' 
#' @export


ven_create = function(nListings, vecCounts, listNames, 
                      locName, groupColors) {
  
  #function can handle only 2 or 3 listings
  if (nListings == 2) {
    
    #coordinates for the circles that compose the venn diagram
    df.venn = data.frame(x = c(0.866, -0.866),
                         y = c(1, 1),
                         labels = c(listNames[1], listNames[2]))
    
    #coordinates for the text that shows the listing counts
    x = c(1.4,-1.4, 0, 2.2, -2.2)
    y = c(1.2,1.2,1.2, 2.1, 2.1)
    
    df.counts = data.frame(x, y, vecCounts)
    
    
              #plot the circles, fill them based on the labels
              #which are tied to a designated color
    ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
      geom_circle(alpha = .3, size = 0.8, colour = 'grey') +
      coord_fixed() +
      theme(legend.title = element_blank()) +
      #group colors designates the fill for each circle
      scale_fill_manual(values = groupColors, name = "") +
      theme_void() + 
      #title is the location
      labs(title = locName) + 
      #place the text
      annotate("text", x = df.counts$x, y = df.counts$y,
               label = df.counts$vecCounts, size = 5)
    
    
  }
  
  else if (nListings == 3) {
    
    #coordinates for the circles
    df.venn <- data.frame(x = c(0, 0.866, -0.866),
                          y = c(1, -0.5, -0.5),
                          labels = c(listNames[1], listNames[2], listNames[3]))
    
    #coordinates for the text that shows the counts for each listing
    x = c(0.05,1.3, -1.3,0.8,-0.8, 0, 0, -0.8, 2, -2)
    y = c(1.5,-0.8, -0.8, 0.6, 0.6, 0.15, -0.8, 2.5,  -1.7, -1.7)
    df.counts = data.frame(x, y, vecCounts)
    
    ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
      geom_circle(alpha = .3, size = 1, colour = 'grey') +
      coord_fixed() +
      scale_fill_manual(values = groupColors, name = "") +
      theme(legend.title = element_blank()) +
      theme_void() + 
      labs(title = locName) +
      annotate("text", x = df.counts$x, y = df.counts$y,
               label = df.counts$vecCounts, size = 5)
    
    
  }
  
  else {
    print("At the moment this function only creates diagrams for 2 or 3 listings")
  }
  
}