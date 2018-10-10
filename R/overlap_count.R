#' Get the counts of individuals who participated in multiple listings, when known.
#' 
#' @description Take overlapdata.csv file and get counts for 
#' @description individuals who participated in multiple listings.
#' 
#' @param overlapData csv file in the format of overlapdata.csv.
#' @param survLoc character object with the corresponding location name in overlapdata.csv.
#' @param listing character object with the corresponding listing name in overlapdata.csv
#' 
#' @return a numeric object
#' 
#' @import dplyr
#' @importFrom magrittr %>%
#' 
#' @export


#in this specific case the survey is the most complete source of data, 
#and this function gets the overlap between the survey one of three
#other types of listings: a  uid, coupon, or rnb
overlap_count = function(survLoc, listing, overlapData) {
  
  if (listing == "uid") {
    
    overlapCount = overlapData %>%
                filter(location == survLoc, uid == "Yes") %>%
                                summarise(n = sum(n)) %>%
                                          pull
  }
  
  else if (listing == "coupon") {
    
    overlapCount = overlapData %>% 
                filter(location == survLoc, coupon == "Yes") %>%
                                summarise(n = sum(n)) %>%
                                          pull
  }
  
  else {
    
    overlapCount = overlapData %>% 
                filter(location == survLoc, rainbow == "Yes") %>%
                                summarise(n = sum(n)) %>%
                                          pull
  }
}