% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/location_sim.R
\name{location_sim}
\alias{location_sim}
\title{Simulates draws from multinomial model of listings for 2 or 3 listings}
\usage{
location_sim(nListings, size, probs)
}
\arguments{
\item{nListings}{number of listings, 2 or 3}

\item{size}{total size drawn}

\item{probs}{vector of inclusion probabilities for each listing}
}
\value{
a list of draws
}
\description{
works for 2 or 3 listings, probs are the inclusion probabilites

for each listing
}
