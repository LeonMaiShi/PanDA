\name{g_FDR}
\alias{g_FDR}
\title{
  Calculating False Discovery Rates (FDR)
}
\description{
\verb{g_FDR} converts the local false discovery rates (fdr) (Efron, 2004) into global FDR by following Newton \emph{et al.} (2004).
}
\usage{
g_FDR(fdr)
}
\arguments{
  \item{fdr}{
     numeric matrix, the local false discovery rates obtained from \code{\link[PanDM]{PanDM}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The global FDR matrix.
}
\references{
Efron, B. (2004). Large-scale simultaneous hypothesis testing: The choice of a null hypothesis. \emph{J. Am. Statist. Assoc.}, \strong{99}(465), 96–104. DOI: \href{http://dx.doi.org/10.1198/016214504000000089}{10.1198/016214504000000089}

Newton, M. A., Noueiry, A., Sarkar, D., and Ahlquist, P. (2004). Detecting differential gene expression with a semiparametric hierarchical mixture method. \emph{Biostatistics}, \strong{5}(2), 155–176. DOI: \href{http://dx.doi.org/10.1093/biostatistics/5.2.155}{10.1093/biostatistics/5.2.155}

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
fdr <- matrix(rnorm(100),10)
FDR <- g_FDR(fdr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
