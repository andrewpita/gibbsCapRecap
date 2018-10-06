
library(dplyr)


#function for extracting the number of individuals who participated 
#in a given listing in a specific location. Assumes input
#survData is in the same format as the regionaldata.csv file

loc_listing_count = function(survLoc, listing, survData) {
  #survLoc is the location name, listing is the specific name
  #of the listing you want a count of, and survData is the 
  #regionaldata.csv file that contains this information
  
  listingCount = survData %>% 
                      filter(location == survLoc) %>%
                              select(listing) %>%
                                    pull
  
}