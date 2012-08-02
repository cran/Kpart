\name{lm.Kpart}
\alias{lm.Kpart}
\title{
something
}
\description{
Offers overall cubic spline model model summary.
}
\usage{
lm.Kpart(data, K)
}
\arguments{
  \item{data}{ Two column matrix, the first column is respone, the second column is predictor in which to base knots.   }
  \item{K}{  A positive interger that will partition the data K times.  }
}
\details{
For each kth partition potential knot selected, linear, quadratic, and cubic terms, use of the regsubsets() function from package(leaps) is implemented to determine the model whch has the lowest BIC. Then after identifying which predictors these are, fits an overall model using lm().
}
\value{
A lm() summary.
}
\references{
Tomas Fitzgerald
}
\author{
Eric Golinko
}
\note{
thank you.
}

\examples{
# example()
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manpi }

