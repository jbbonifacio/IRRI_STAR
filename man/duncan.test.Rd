\name{duncan.test}
\alias{duncan.test}
%- duncan.
\title{ Duncan's new multiple range test }
\description{
  This test is adapted from the Newman-Keuls method. Duncan's test does not control family wise error rate at the specified alpha level. It has more power than the other post tests, but only because it doesn't control the error rate properly. The Experimentwise Error Rate at: 1-(1-alpha)^(a-1); where "a" is the number of means and is the Per-Comparison Error Rate. Duncan's procedure is only very slightly more conservative than LSD. The level by alpha default is 0.05.
}
\usage{
duncan.test(y, trt, DFerror, MSerror, alpha = 0.05, group=TRUE, main = NULL, pwOrder = c("trmt","means"))
}

\arguments{
  \item{y}{ model(aov or lm) or answer of the experimental unit  }
  \item{trt}{ Constant( only y=model) or vector treatment applied to each experimental unit  }
  \item{DFerror}{ Degree free }
  \item{MSerror}{ Mean Square Error }
  \item{alpha}{ Significant level }
  \item{group}{ TRUE or FALSE }
  \item{main}{ Title }
  \item{pwOrder}{means or trmt}
}
\details{
  It is necessary first makes a analysis of variance.
}
\value{

  \item{y }{class (aov or lm) or vector numeric}
  \item{trt }{constant (only y=model) or vector alfanumeric}
  \item{DFerror }{Numeric}
  \item{MSerror }{Numeric}
  \item{alpha }{Numeric}
  \item{group}{ Logic}  
  \item{main }{Text}

}
\references{ Principles and procedures of statistics a biometrical approach
             Steel & Torry & Dickey. Third Edition 1997 }
\author{ Felipe de Mendiburu of package agricolae}

\seealso{\code{\link{LSD.test}}, \code{\link{HSD.test} }, \code{\link{SNK.test} } }

\keyword{ htest }

