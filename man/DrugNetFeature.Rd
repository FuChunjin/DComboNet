% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DrugNetFeatures.R
\name{DrugNetFeature}
\alias{DrugNetFeature}
\title{Drug drug interaction network edge weight: Drug-Drug similarity}
\usage{
DrugNetFeature(
  druglist = NULL,
  cellline = NULL,
  model = c("L1", "L2"),
  CDK_FP = NULL,
  pubchemFP = NULL,
  MACCS_FP = NULL,
  load_dir
)
}
\arguments{
\item{druglist}{drug and its corresponding "drugbank ID" in \code{data
frame}}

\item{cellline}{cancer cell lines name, provide only when level two model is
called (\code{model = "L2"}).}

\item{model}{define which level model should be provided. Options are level
one model ("\code{L1}") and level two model ("\code{L2}")}

\item{CDK_FP}{Drug fingerprint matrix generated from PEDAL}

\item{pubchemFP}{Drug Pubchem fingerprint matrix generated from PEDAL}

\item{MACCS_FP}{Drug MACCS fingerprint matrix generated from PEDAL}

\item{load_dir}{path to load and save modeling data files}
}
\value{
\code{features} table that includes the three drug similarity values
  and the integrated pharmacological score
}
\description{
The function \code{DrugNetFeature} is to calculate drug ATC
  code-based similarity, chemical structure fingerprint-based similarity and
  drug side effect-based similarity and integrate together as the final
  pharmacological similarity.
}
\details{
\code{DrugNetFeature} integrates calculation of three drug-drug
  similarity functions, \code{druglist} included "drugbank ID" is
  required to provide for file preparation. drug ATC code-based similarity,
  chemical structure fingerprint-based similarity and drug side effect-based
  similarity will be computed seperately.

\code{NA} generated due to lack of data will be imputed by the \code{mean}
  of each similarity score. Three simiarity score will be integrated into
  phamacological score and the score of known combinations will be replaced
  as \code{1}. Final \code{feature} table will be saved in provided
  \code{load_dir} path under \emph{/drug_net/} folder named with
  "features.csv" for level one model or "features_\emph{cellline}.csv" if
  level two model is called and \code{cellline} is provided.
}
\examples{

load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
druglist = NULL
# or druglist = read.table("druglist_example.csv", sep = ",", header = TRUE,
# stringsAsFactors = FALSE)
# For level one model
cellline = NULL
model = "L1"

\dontrun{
# For level two model
cellline = "HEPG2"
model = "L2"
}

# If no newly inputed drug(s): druglist = NULL
CDK_FP=NULL
pubchemFP=NULL
MACCS_FP=NULL
# Otherwise, please input PaDEL generated three fingerprint files
\dontrun{DrugNetFeature(druglist = druglist, cellline = cellline, model = model,
CDK_FP = CDK_FP, pubchemFP = pubchemFP, MACCS_FP = MACCS_FP, load_dir =
load_dir)}

}
\seealso{
\code{\link{ATC.sim}} for ATC code based similarity,
  \code{\link{structure.sim}} for chemical fingerprint based similarity,
  \code{\link{SE.sim}} for side effect based similarity.
}
