% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf2weight.R
\name{rf2weight}
\alias{rf2weight}
\title{Predict a candidate gene for each region conditioning on the candidates selected in rest regions.}
\usage{
rf2weight(pickup, region_chosen, remaining, opt, count)
}
\arguments{
\item{pickup}{Numeric. Region index.
region_chosen A vector of TRUE/FALSE indicates which regions are selected as conditioning regions
remaining A vector of genes chosen in the last iteration ((n-1)th round).
opt A list of parameter in Gibbs}
}
\value{
Characters of a single gene
}
\description{
Predict a candidate gene for each region conditioning on the candidates selected in rest regions.
}
\examples{
predict_gene(pickup,region_chosen,remaining,opt,count))

}
