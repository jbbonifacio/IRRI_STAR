\name{cvm.test}
\alias{cvm.test}
\title{Cramer-von Mises test for normality}
\description{
  Performs the Cramer-von Mises test for the composite hypothesis of normality, 
  see e.g. Thode (2002, Sec. 5.1.3).
}

\usage{
cvm.test(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data values, the number of 
  which must be greater than 7. Missing values are allowed.}
}
\details{The Cramer-von Mises test is an EDF omnibus test for the composite hypothesis of normality. 
The test statistic is 
\deqn{
W = \frac{1}{12 n} + \sum_{i=1}^{n} (p_{(i)} - \frac{2i-1}{2n}),
}
where \eqn{p_{(i)} = \Phi([x_{(i)} - \overline{x}]/s)}. Here, 
\eqn{\Phi} is the cumulative distribution function 
of the standard normal distribution, and \eqn{\overline{x}} and \eqn{s} 
are mean and standard deviation of the data values. 
The p-value is computed from the modified statistic 
\eqn{Z=W (1.0 + 0.5/n)} according to Table 4.9 in 
Stephens (1986). 
}
\value{
  A list with class "htest" containing the following components: 
  \item{statistic}{the value of the Cramer-von Mises statistic.}
  \item{p.value }{the p-value for the test.}
  \item{method}{the character string "Cramer-von Mises normality test".}
  \item{data.name}{a character string giving the name(s) of the data.} 
}
\references{Stephens, M.A. (1986): Tests based on EDF statistics. In:
D'Agostino, R.B. and Stephens, M.A., eds.: Goodness-of-Fit Techniques.
Marcel Dekker, New York.

Thode Jr., H.C. (2002): Testing for  Normality. Marcel Dekker, New York.
}
\author{Juergen Gross from the package "nortest"}


\seealso{\code{\link{shapiro.test}} for performing the Shapiro-Wilk test for normality. 
\code{\link{ad.test}}, \code{\link{lillie.test}}, 
\code{\link{sf.test}} for performing further tests for normality.
\code{\link{qqnorm}} for producing a normal quantile-quantile plot.}

\examples{
cvm.test(rnorm(100, mean = 5, sd = 3))
cvm.test(runif(100, min = 2, max = 4))

}
\keyword{htest}
