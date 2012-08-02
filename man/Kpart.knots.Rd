\name{Kpart.knots}
\alias{Kpart.knots}
\title{
Returns the spline knots selected by lm.Kpart. 
}
\description{
Displays selected knots by best subsets method after each absolute maximum deviate in partition mean is determined.
}
\usage{
Kpart.knots(data, K)
}
\arguments{
  \item{data}{ Two column matrix, the first column is respone, the second column is predictor in which to base knots.  }
  \item{K}{ A positive interger that will partition the data K times. }
}
\details{
Use of absolute maximum deviate from partition mean. 
}
\value{
Will return the values of the second column(predictor) in which the knots were selected from lm.Kpart.
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
\keyword{ manip }
