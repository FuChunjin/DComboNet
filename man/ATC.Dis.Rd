% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DrugNetFeatures.R
\name{ATC.Dis}
\alias{ATC.Dis}
\title{Drug ATC Code Distance Matrix Function}
\usage{
ATC.Dis(drug_atc, load_dir)
}
\arguments{
\item{drug_atc}{a data frame of drug and its corresponding ATC code
in ATC system, generated from function \code{\link{drugATCprep}}}

\item{load_dir}{path to load or save modelling data files}
}
\value{
\code{ATC.Dis} returns ATC code distance matrix, an intermediate file
  for computing drug-drug ATC code based similarity
}
\description{
The function \code{ATC.Dis} is to calculate drug ATC code distance in
  Anatomical Therapeutic Chemical Classification System.
}
\examples{

druglist = data.frame(drug = "Sorafenib")
load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
drug_atc = drugATCprep(druglist = druglist, load_dir = load_dir)
ATC_dis = ATC.Dis(drug_atc = drug_atc, load_dir = load_dir)

}
\seealso{
\code{\link{drugATCprep}} to prepare input drug-atc file
}
