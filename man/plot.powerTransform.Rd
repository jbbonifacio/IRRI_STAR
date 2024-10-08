\name{plot.powerTransform}
\Rdversion{1.1}
\alias{plot.powerTransform}
\title{ plot Method for powerTransform Objects
}
\description{
This function provides a simple function for plotting data using power 
transformations.

}
\usage{
\method{plot}{powerTransform}(x, z = NULL, round = TRUE, plot = pairs, ...)
}

\arguments{
  \item{x}{name of the power transformation object
}
  \item{z}{
  Additional variables of the same length as those used to get the transformation
  to be plotted, default is \code{NULL}.
}
  \item{round}{
  If \code{TRUE}, the default, use rounded transforms, if \code{FALSE} use the MLEs.
}
  \item{plot}{
  Plotting method.  Default is \code{pairs}.  Another possible choice is
  \code{scatterplot.matrix} from the \code{car} package.
}
  \item{\dots}{
  Optional arguments passed to the plotting method
}
}
\details{
The data used to estimate transformations using \code{powerTransform} are plotted
in the transformed scale.
}
\value{
None.  Produces a graph as a side-effect.
}
\references{
Weisberg, S. (2005) \emph{Applied Linear Regression}, Third Edition. Wiley.

Fox, J. and Weisberg, S. (2011) \emph{An R Companion to Applied Linear Regression},
Second Edition, Sage.
}
\author{ Sanford Weisberg, <sandy@stat.umn.edu> }


\seealso{
\code{\link{powerTransform}}
}
\keyword{ regression }% __ONLY ONE__ keyword per line
