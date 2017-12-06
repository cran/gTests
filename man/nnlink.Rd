\name{nnlink}
\alias{nnlink}
\title{Construct similarity graph by 1-NNL}
\description{This function provides the edges of the similarity graph constructed by 1-NNL.}
\usage{
nnlink(distance)
}
\arguments{
  \item{distance}{A K by K matrix, which is the distance matrix on the distinct values and K is the number of distinct values with at least one observation in either group.}
}
\value{
  \item{E}{An edge matrix representing a similarity graph on the distinct values with the number of edges in the similarity graph being the number of rows and 2 columns.  Each row records the subject indices of the two ends of an edge in the similarity graph.}
}
\seealso{
  \code{\link{getGraph}}
}

