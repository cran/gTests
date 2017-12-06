\name{g.tests}
\alias{g.tests}
\title{Graph-based two-sample tests}
\description{This function provides four graph-based two-sample tests.}
\usage{
g.tests(E, sample1ID, sample2ID, test.type="all", maxtype.kappa = 1.14, perm=0)
}
\arguments{
  \item{E}{An edge matrix representing a similarity graph with the number of edges in the similarity graph being the number of rows and 2 columns. Each row records the subject indices of the two ends of an edge in the similarity graph.}
  \item{sample1ID}{The subject indices of sample 1.}
  \item{sample2ID}{The subject indices of sample 2.}
  \item{test.type}{The default value is "all", which means all four tests are performed: orignial edge-count test (Friedman and Rafsky (1979)), generalized edge-count test (Chen and Friedman (2016)), weighted edge-count test (Chen, Chen and Su (2016)) and maxtype edge-count tests (Zhang and Chen (2017)). Set this value to "original" or "o" to permform only the original edge-count test; set this value to "generalized" or "g" to perform only the generalized edge-count test; set this value to "weighted" or "w" to perform only the weighted edge-count test; and set this value to "maxtype" or "m" to perform only the maxtype edge-count tests.}
  \item{maxtype.kappa}{The value of parameter(kappa) in the maxtype edge-count tests. The default value is 1.14.}
  \item{perm}{The number of permutations performed  to calculate the p-value of the test. The default value is 0, which means the permutation is not performed and only approximate p-value based on asymptotic theory is provided. Doing permutation could be time consuming, so be cautious if you want to set this value to be larger than 10,000.}
}
\value{
  \item{test.statistic}{The test statistic.}
  \item{pval.approx}{The approximated p-value based on asymptotic theory.}
  \item{pval.perm}{The permutation p-value when argument `perm` is positive.}
}
\references{
Friedman J. and Rafsky L. Multivariate generalizations of the WaldWolfowitz
and Smirnov two-sample tests. The Annals of Statistics, 7(4):697-717, 1979.

Chen, H. and Friedman, J. H.  A new graph-based two-sample test for multivariate and object data. Journal of the American Statistical Association, 2016.

Chen, H., Chen, X. and Su, Y. A weighted edge-count two sample test for multivariate and object data. Journal of the American Statistical Association, 2017.

Zhang, J. and Chen, H. Graph-based two-sample tests for discrete data.
}
\examples{
# the "example" data contains three similarity graphs represted in the matrix form: E1, E2, E3.
data(example) 
 
# E1 is an edge matrix representing a similarity graph.
# It is constructed on two samples with mean difference.
# Sample 1 indices: 1:100; sample 2 indices: 101:250.
g.tests(E1, 1:100, 101:250) 

# E2 is an edge matrix representing a similarity graph.
# It is constructed on two samples with variance difference.
# Sample 1 indices: 1:100; sample 2 indices: 101:250.
g.tests(E2, 1:100, 101:250)

# E3 is an edge matrix representing a similarity graph.
# It is constructed on two samples with mean and variance difference.
# Sample 1 indices: 1:100; sample 2 indices: 101:250.
g.tests(E3, 1:100, 101:250)

## Uncomment the following line to get permutation p-value with 200 permutations.
# g.tests(E1, 1:100, 101:250, perm=200)
}

