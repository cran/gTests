\name{g.tests_discrete}
\alias{g.tests_discrete}
\title{Graph-based two-sample tests for discrete data}
\description{This function provides four graph-based two-sample tests for discrete data.}
\usage{
g.tests_discrete(E, counts, test.type = "all", maxtype.kappa = 1.14, perm = 0)
}
\arguments{
  \item{E}{An edge matrix representing a similarity graph on the distinct values with the number of edges in the similarity graph being the number of rows and 2 columns. Each row records the subject indices of the two ends of an edge in the similarity graph.}
  \item{counts}{A K by 2 matrix, where K is the number of distinct values. It specifies the counts in the K distinct values for the two samples.}
  \item{test.type}{The default value is "all", which means all four tests are performed: the orignial edge-count test (Chen and Zhang (2013)), extension of the generalized edge-count test (Chen and Friedman (2016)), extension of the weighted edge-count test (Chen, Chen and Su (2016)) and extension of the maxtype edge-count tests (Zhang and Chen (2017)).  Set this value to "original" or "o" to permform only the original edge-count test; set this value to "generalized" or "g" to perform only extension of the generalized edge-count test; set this value to "weighted" or "w" to perform only extension of the weighted edge-count test; and set this value to "maxtype" or "m" to perform only extension of the maxtype edge-count tests.}
  \item{maxtype.kappa}{The value of parameter(kappa) in the extension of the maxtype edge-count tests. The default value is 1.14.}
  \item{perm}{The number of permutations performed  to calculate the p-value of the test. The default value is 0, which means the permutation is not performed and only approximate p-value based on asymptotic theory is provided. Doing permutation could be time consuming, so be cautious if you want to set this value to be larger than 10,000.}
}
\value{
  \item{test.statistic_a}{The test statistic using `average` method to construct the graph.}
  \item{test.statistic_u}{The test statistic using `union` method to construct the graph.}
  \item{pval.approx_a}{Using `average` method to construct the graph, the approximated p-value based on asymptotic theory.}
  \item{pval.approx_u}{Using `union` method to construct the graph, the approximated p-value based on asymptotic theory.}
  \item{pval.perm_a}{Using `average` method to construct the graph, the permutation p-value when argument `perm` is positive.}
  \item{pval.perm_u}{Using `union` method to construct the graph, the permutation p-value when argument `perm` is positive.}
}
\references{
Friedman J. and Rafsky L. Multivariate generalizations of the WaldWolfowitz
and Smirnov two-sample tests. The Annals of Statistics, 7(4):697-717, 1979.

Chen, H. and Zhang, N. R. Graph-based tests for two-sample comparisons of categorical data. Statistica Sinica, 2013.

Chen, H. and Friedman, J. H.  A new graph-based two-sample test for multivariate and object data. Journal of the American Statistical Association, 2016.

Chen, H., Chen, X. and Su, Y. A weighted edge-count two sample test for multivariate and object data. Journal of the American Statistical Association, 2017.

Zhang, J. and Chen, H. Graph-based two-sample tests for discrete data.
}
\examples{
# the "example_discrete" data contains three two-sample counts data 
# represted in the matrix form: counts1, counts2, counts3 
# and the corresponding distance matrix on the distinct values: ds1, ds2, ds3.
data(example_discrete) 

# counts1 is a K by 2 matrix, where K is the number of distinct values. 
# It specifies the counts in the K distinct values for the two samples. 
# ds1 is the corresponding distance matrix on the distinct values. 
# The data is generated from two samples with mean shift.
Knnl = 3
E1 = getGraph(counts1, ds1, Knnl, graph = "nnlink")
g.tests_discrete(E1, counts1)
 
# counts2 is a K by 2 matrix, where K is the number of distinct values. 
# It specifies the counts in the K distinct values for the two samples. 
# ds2 is the corresponding distance matrix on the distinct values. 
# The data is generated from two samples with spread difference.
Kmst = 6
E2 = getGraph(counts2, ds2, Kmst, graph = "mstree")
g.tests_discrete(E2, counts2)
 
# counts3 is a K by 2 matrix, where K is the number of distinct values. 
# It specifies the counts in the K distinct values for the two samples. 
# ds3 is the corresponding distance matrix on the distinct values. 
# The data is generated from two samples with mean shift and spread difference.
Knnl = 3
E3 = getGraph(counts3, ds3, Knnl, graph = "nnlink")
g.tests_discrete(E3, counts3)

## Uncomment the following line to get permutation p-value with 200 permutations.
# Knnl = 3
# E1 = getGraph(counts1, ds1, Knnl, graph = "nnlink")
# g.tests_discrete(E1, counts1, test.type = "all", maxtype.kappa = 1.31, perm = 300)
}








