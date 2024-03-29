% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PGE_model.R
\name{nls_3peak_two_tissue_PGE}
\alias{nls_3peak_two_tissue_PGE}
\title{Fitting two indipendent peaks to the kmer histogram, this model aims to allow for testing of evenness of spacing (expected 1:2)}
\usage{
nls_3peak_two_tissue_PGE(
  x,
  y,
  k,
  estKmercov,
  estLength,
  estR,
  biasEst,
  estFrac,
  max_iterations
)
}
\arguments{
\item{x, }{y vectors with coverages and theuir respective numbers of kmers (frequencies)}

\item{estKmercov, }{estLength, biasEst parameters that serve as initial estimates of kmer coverage, genome length and overdispersal}

\item{estR}{- heterozygosity, a fixed parameter treated as a constant}
}
\description{
Fitting two indipendent peaks to the kmer histogram, this model aims to allow for testing of evenness of spacing (expected 1:2)
}
