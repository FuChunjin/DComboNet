% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RWR_fun.R
\name{rwr}
\alias{rwr}
\title{Random walk with restart algorithm}
\usage{
rwr(tm, r = 0.7, seeds_score)
}
\arguments{
\item{tm}{transition matrix generated from drug-gene/pathway interaction
networks via \code{TransitionMatrix}}

\item{r}{numeric, global restart parameter}

\item{seeds_score}{drugseed-score \code{data.frame}, generated via
\code{seedscore}}
}
\value{
global proximity between drug seed and other drug in
  drug-gene/pathway interaction network.
}
\description{
The function \code{rwr} is to calculate global proximity between
  drug seed and other drugs in complex network.
}
\details{
The equation of random walk with restart is as follow:
  \deqn{P_{t+1} = (1-\sigma)H * P_{t} + \sigma P_{0}}
}
\examples{
\dontrun{
tm = TransitionMatrix(...)
seeds = seedscore(seeds = "Sorafenib")
rwr(tm, r = 0.7, seeds_score = seeds)
}

}
\references{
Valdeolivas A, Tichit L, Navarro C, et al. Random walk with
  restart on multiplex and heterogeneous biological networks[J].
  Bioinformatics, 2019, 35(3): 497-505.
}
