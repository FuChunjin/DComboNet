% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenePathwayNet.R
\name{genepathwayAdj}
\alias{genepathwayAdj}
\title{Generate the bipartite matrix from gene-pathway interaction networks}
\usage{
genepathwayAdj(drugAdj, pathwayadj, geneadj, gpweight, load_dir)
}
\arguments{
\item{drugAdj}{drug adjacency matrix (sparse matrix) generated from drug-drug
association network via \code{AdjMatrix_drugs}}

\item{pathwayadj}{pathway-pathway association adjacency matrix (sparse
matrix)}

\item{geneadj}{gene-gene association network adjacency matrix (sparse matrix)}

\item{gpweight}{the weight of gene-pathway edges}

\item{load_dir}{path to load or save modeling data files}
}
\value{
\code{genepathwayAdj} gene-pathway adjacency matrix
}
\description{
The function \code{genepathwayAdj} is to generate the adjacency
  matrix from gene-pathway interaction network
}
\examples{

}
