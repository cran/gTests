\name{nnlink_K}
\alias{nnlink_K}
\title{Construct similarity graph by k-NNL}
\description{This function provides the edges of the similarity graph constructed by k-NNL.}
\usage{
nnlink_K(distance,K)
}
\arguments{
  \item{distance}{A K by K matrix, which is the distance matrix on the distinct values and K is the number of distinct values with at least one observation in either group.}
  \item{K}{Set the value of k in "k-NNL" to construct the similarity graph.}
}
\value{
  \item{E}{An edge matrix representing a similarity graph on the distinct values with the number of edges in the similarity graph being the number of rows and 2 columns.  Each row records the subject indices of the two ends of an edge in the similarity graph.}
}
\seealso{
  \code{\link{getGraph}}
}

