% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_deg_groups.R
\name{get_deg_groups}
\alias{get_deg_groups}
\title{Create FC groups with counts based on the number of genes wanted}
\usage{
get_deg_groups(n_genes, fc_range_min, fc_range_max, seed, deg_fc_interval)
}
\arguments{
\item{fc_range_min}{a number between 1 and 10. The minimum amount of fold 
change to simulate a difference in expression between case and control 
groups. This is used to derive the differentially expressed genes. The
fold change changes genes in either direction, so a fold change of 2 would
cause the mean expression in cases to be either twice the amount or half the
amount for any particular gene. Defaults 1.1 which arguably would be hard to
by a model since the FC difference between case and controls is so low..}

\item{fc_range_max}{a number between 1 and 10. The maximum amount of fold 
change to simulate a difference in expression between case and control 
groups. This is used to derive the differentially expressed genes. The
fold change changes genes in either direction, so a fold change of 2 would
cause the mean expression in cases to be either twice the amount or half the
amount for any particular gene. Defaults 3. n_genes will be randomly 
sampled between fc_range_min and fc_range_max at an interval of 0.1.}

\item{seed}{Use this to ensure the same count datasets are generated in 
subsequent runs, this is helpful when comparing models so we know the exact
same dataset was used. Default False, so no seed is set. If integer value
passed, seed is set using it.}

\item{deg_fc_interval}{the incremental value between fc_range_min and 
fc_range_max used to derive the number of fold change groups. This uniform 
sampling of fold change values is used to create DEGs from a distribution 
with the specified fold change. Default is 0.5}
}
\value{
named vector containing the FC value and the number of genes to be 
derived for it
}
\description{
Create FC groups with counts based on the number of genes wanted
}
\keyword{internal}
