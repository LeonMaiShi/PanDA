\name{PanDM_BIC}
\alias{PanDM_BIC}
\title{BIC for PanDM
}
\description{
Generate the Bayesian Information Criterion (BIC) table and plot after \code{\link[PanDM]{PanDM}} model fitting.
}
\usage{
PanDM_BIC(resultlist, BICplot = TRUE)
}
\arguments{
  \item{resultlist}{
        a list of model fitting results from \emph{PanDM} for different \emph{K}s.
}
  \item{BICplot}{
        logical; if TRUE, then the relationship between BIC and the number of clusters (\emph{K}) will be plotted.
}
}
\details{
The Bayesian Information Criterion (BIC) is a criterion for model selection among a finite set of distinct models (Schwarz, 1978). We treat the \emph{PanDM} fitting for each \emph{K} as a different model from others. The model with the lowest BIC is adopted as the best one. The corresponding \emph{K} is regarded as the optimal number of DM patterns in the input dataset of \emph{PanDM}.
}
\value{
The BIC table, as well as BIC plot if \verb{BICplot = TRUE}.
}
\references{
Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}, \strong{6}(2), 461-464. DOI: \href{http://projecteuclid.org/euclid.aos/1176344136}{10.1214/aos/1176344136}

Shi, M., Tsui, S., Wu, H., Wei, Y.,(2017) \emph{Pan-cancer analysis of differential DNA methylation patterns}
}
\author{
Mai Shi, Stephen Kwok-Wing Tsui, Hao Wu, Yingying Wei
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{PanDM_BIC(resultlist, BICplot = TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
