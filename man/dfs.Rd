\name{dfs}
\alias{dfs}
\title{Depth-first search}
\description{ One starts at the root and explores as far as possible along each branch before backtracking.}
\usage{
dfs(s,visited,adj)
}
\arguments{
  \item{s}{The root node.}
  \item{visited}{N by 1 vector, where N is the number of nodes. This vector records whether nodes have been visited or not with 1 if visited and 0 otherwise.}
  \item{adj}{N by N adjacent matrix.}
}
\seealso{
  \code{\link{getGraph}}
}

