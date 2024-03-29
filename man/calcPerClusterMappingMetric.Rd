% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confidenceScores.R
\name{calcPerClusterMappingMetric}
\alias{calcPerClusterMappingMetric}
\title{Per-cluster Confidence Score:
Calculates the Mahalanobis distance from user-defined query clusters to their nearest
reference centroid after initial projection into reference PCA space.
All query cells in a cluster get the same score. Higher distance indicates less confidence.
Due to the instability of estimating covariance with small numbers of cells, we do not assign a
score to clusters smaller than u * d, where d is the dimensionality of the embedding and u is specified.}
\usage{
calcPerClusterMappingMetric(
  reference,
  query,
  query_cluster_labels,
  metric = "mahalanobis",
  u = 2,
  lambda = 0
)
}
\arguments{
\item{reference}{Reference object as returned by Symphony buildReference()}

\item{query}{Query object as returned by Symphony mapQuery()}

\item{query_cluster_labels}{Vector of user-defined labels denoting clusters / putative novel cell type to calculate the score for}

\item{metric}{Uses Mahalanobis by default, but added as a parameter for potential future use}

\item{u}{Do not assign scores to clusters smaller than u * d (see above description)}

\item{lambda}{Optional ridge parameter added to covariance diagonal to help stabilize numeric estimates}
}
\value{
A data.frame of per-cluster mapping metric scores for each user-specified query cluster.
}
\description{
Per-cluster Confidence Score:
Calculates the Mahalanobis distance from user-defined query clusters to their nearest
reference centroid after initial projection into reference PCA space.
All query cells in a cluster get the same score. Higher distance indicates less confidence.
Due to the instability of estimating covariance with small numbers of cells, we do not assign a
score to clusters smaller than u * d, where d is the dimensionality of the embedding and u is specified.
}
