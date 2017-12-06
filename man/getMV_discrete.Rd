\name{getMV_discrete}
\alias{getMV_discrete}
\title{Get intermediate results for g.tests_discrete function}
\description{This function calculates means and variances of R1 and R2 quantities using `average` method and `union` method to construct the graph.}
\usage{
getMV_discrete(E,vmat)
}
\arguments{
  \item{E}{An edge matrix representing a similarity graph on the distinct values with the number of edges in the similarity graph being the number of rows and 2 columns.  Each row records the subject indices of the two ends of an edge in the similarity graph.}
  \item{vmat}{A K by 2 matrix, where K is the number of distinct values with at least one observation in either group. It specifies the counts in the K distinct values for the two samples.}
}
\seealso{
  \code{\link{g.tests_discrete}}
}

