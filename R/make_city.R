#' Get the values for parameters for each location for use in simulation.
#' 
#' @description To simulate data we need values for parameters at each location 
#' 
#' @param string The output of the Gibbs Chain contains parameter estimates labeled by location. This string corresponds to that location. 
#' @param df a named vector or df that contains the parameter values
#' 
#' @return a numeric vector
#' 


make_city = function(string, df) {
  
  indeces = grepl(paste(string,"$", sep = ""), names(df))
  
  city = df[indeces]
  names(city) = NULL
  return(city)
  
}

