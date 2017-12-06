\name{getGraph}
\alias{getGraph}
\title{Construct similarity graph}
\description{This function provides two methods to construct the similarity graph.}
\usage{
getGraph(counts, mydist, K, graph.type = "mstree")
}
\arguments{
  \item{counts}{A K by 2 matrix, where K is the number of distinct values. It specifies the counts in the K distinct values for the two samples.}
  \item{mydist}{A K by K matrix, which is the distance matrix on the distinct values.}
  \item{K}{Set the value of k in "k-MST" or "k-NNL" to construct the similarity graph.}
  \item{graph.type}{Specify the type of the constructing graph. The default value is "mstree", which means constructing the minimal spanning tree as the similarity graph. Set this value to "nnlink" to construct the similarity graph by the nearest neighbor link method.}
}
\value{
  \item{E}{An edge matrix representing a similarity graph on the distinct values with the number of edges in the similarity graph being the number of rows and 2 columns.  Each row records the subject indices of the two ends of an edge in the similarity graph.}
}
\seealso{
  \code{\link{g.tests_discrete}}
}

