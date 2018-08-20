\name{PanDM}
\alias{PanDM}
\title{Pan-cancer Differential Methylation analysis
}
\description{
This function fits a hierarchical two-normal-component mixture model to the summary statistics resulted from differential methylation analyses for multiple cancer types. It uses expectation-maximization (EM) algorithm to estimate parameters and characterizes the data with a differential methylation pattern matrix and a proportion vector corresponding to each cluster in the matrix. The probabilities of being differentially methylated in every cancer type is also given in the results.
}
\usage{
PanDM(dat, K, epsilon = 1e-04, max.iter = 5000, save.initial = TRUE, para.cores = 0)
}
\arguments{
  \item{dat}{
        input data matrix composed of the \emph{z-values} transformed transformed from \emph{p-values}. Each row subject can be a probe (CpG site), a gene or a differentially methylated region (DMR). Each column corresponds to a cancer type.
  }
  \item{K}{
        a scalar, the number of clusters assumed existing in the data.
  }
  \item{epsilon}{
        a scalar, the tolerance bound of error served as convergence criterion. The EM algorithm is considered to reach convergence and the program is subsequently ended when the maximum absolute change of the parameters is less than \verb{epsilon}.
  }
  \item{max.iter}{
        the maximum number of iterations for estimating parameters by the EM algorithm.
  }
  \item{save.initial}{
        logical; if TRUE, then the initial values for parameter estimation will be saved in the current working directory.
  }
  \item{para.cores}{
        a scalar, the number of cores used for parallel computing. Both \verb{0} and \verb{1} stand for not using parallel. Please check \strong{Details} for more information regarding this parameter.
  }
}
\details{
PanDM speeds up computation by using the function \code{\link[parallel]{mclapply}} from the \code{\link[parallel]{parallel}} package. Please kindly note that \code{\link[parallel]{mclapply}} relies on forking and thus is not available on Windows unless \verb{para.cores = 0} or \verb{1}.
}
\value{
PanDM returns a list of items:
  \item{cluster}{
        a vector of integers, the cluster membership for each row subject.
  }
  \item{p.post}{
        a matrix composed of the probability for each row subject to be differentially methylated in every cancer type. The orders of rows and columns remain the same as those in the input \verb{dat}.
  }
  \item{p}{
        a vector suggesting the proportion of \verb{K} clusters.
  }
  \item{q}{
        the differential methylation pattern matrix. Each row corresponds to a cluster and each column corresponds to a cancer type.
  }
  \item{loglike}{
        the log-likelihood of the final model.
  }
  \item{bic}{
        the Bayesian Information Criterion (BIC) values of the final model, used to select the optimal number of K.
  }
  \item{mu}{
        a two-row matrix of the estimated component means of the normal-mixture model. Each row corresponds to one of the two components and each column corresponds to a cancer type.
  }
  \item{sigma}{
        a two-row matrix of the estimated component standard deviations of the normal-mixture model. Each row corresponds to one of the two components and each column corresponds to a cancer type.
  }
}
\references{
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
data(simudat)
fit.results <- PanDM(dat = simudat, K = 9, epsilon = 1e-04, max.iter = 10,
                     save.initial = FALSE, para.cores = 0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
