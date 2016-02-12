\name{Concordance-package}
\alias{Concordance-package}
\alias{Concordance}
\alias{concordance.test}
\docType{package}
\title{
Concordance.test
}
\description{
Determines the non-parametric concordance coeffient between two or more replicate measurements.
}

\usage{
concordance.test(mat, mat2 = NULL, alternative =NULL, exact = FALSE)
}
\author{
R.Kuiper
}

\arguments{
  \item{mat}{
A matrix of measurements containing the subjects in the rows and observers in the columns. Missing values can be set to NA.
}
  \item{mat2}{Optional (See mat1); If given, a test for difference between concordances is performed.
}
  \item{alternative}{
The alternative hypothesis. One of : "less", "greater", "two.sided". Default is 'less' for a concordance test and 'two.sided' for a difference between concordances.
}
  \item{exact}{
	A logical value: TRUE, FALSE. In case of a concordance test, an exact test is performed. In case of testing difference between concordances, it is ignorned.
}
 
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##Define two (random) matrices of 10 observation each having 3 replicates
n<-10;b<-3; 
mat1<-matrix(sample(n*b),ncol=b)
mat2<-matrix(sample(n*b),ncol=b)

#Perform an one sided exact test for concordance between replicates:
concordance.test(mat1,exact=TRUE)

#Perform a two sided test for equal concordances:
concordance.test(mat1,mat2)

}

