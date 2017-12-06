\name{permute_discrete}
\alias{permute_discrete}
\title{Generate a permutation for two discrete data groups}
\description{This function permutes the observations maintaining the two sample sizes unchaged.}
\usage{
permute_discrete(vmat)
}
\arguments{
  \item{vmat}{A K by 2 matrix, where K is the number of distinct values with at least one observation in either group. It specifies the counts in the K distinct values for the two samples.}
}
\seealso{
  \code{\link{g.tests_discrete}}
}

