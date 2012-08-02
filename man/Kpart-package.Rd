\name{Kpart-package}
\alias{Kpart-package}
\alias{Kpart}
\docType{package}
\title{
Kpart
}
\description{
Kpart spline fitting
}
\details{
\tabular{ll}{
Package: \tab Kpart\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-08-02\cr
License: \tab Open Source\cr
}
~~ This package is intended for use with non-lineraly associated data. The function lm.Kpart firsts selects points for cubic spline knots using an alogorithm to find the abolute maximum deviate from the partiton mean, then fits a best fitting model by using the best subset method and minimum BIC. The function knots.Kpart returns the values selcted as knots in the model given by lm.Kpart, and plot.Kpart returns a plot with fitted values and vertical lines representing the spline points. Kpart(data,K) takes two arguments, data is a two column matrix where the first column is the repsonse and the second column is the predictor upon which knots are selected. K is a postive interger that indicates how many equally spaced partitions the user would like to produce.~~
}
\author{
Eric Golinko

Maintainer: egolinko@gmail.com
}

\keyword{ package }
