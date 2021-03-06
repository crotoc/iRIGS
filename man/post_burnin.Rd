% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_burnin.R
\name{post_burnin}
\alias{post_burnin}
\title{Two-step sampling bundle consist of a burnin step and a post burnin step.}
\usage{
post_burnin(opt, evi_mode, with_region, run_mode = "training")
}
\arguments{
\item{opt}{A list of parameters in Gibbs script
mode burnin/post_burnin mode
with_region Logical. Whether estimate region risk.}
}
\value{
A list containing results of sampling
}
\description{
Two-step sampling bundle consist of a burnin step and a post burnin step.
}
\examples{
opt$training_gene<-post_burnin(opt,evi_mode="generic",with_region = FALSE,step="training")

}
