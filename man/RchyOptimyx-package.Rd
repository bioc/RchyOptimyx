\name{RchyOptimyx-package}
\alias{RchyOptimyx-package}
\docType{package}
\title{
Gating Hierarchy Optimization for Flow Cytometry
}
\description{
RchyOptimyx is a computational tool which determines the minimal sets of
markers that can identify a target population within a given purity
standard. This can enable panel redesign using smaller subsets of
markers, with consequent savings on reagents, potential for adaptation to older instruments, or the addition of other markers of interest to the panel design.
}
\details{
\tabular{ll}{
Package: \tab RchyOptimyx\cr
Type: \tab Package\cr
Version: \tab 0.99.1\cr
Date: \tab 2011-10-12\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}
See the \code{RchyOptimyx} function.
}

\keyword{ package }
\seealso{
  RchyOptimyx
}

\references{
D. Eppstein, Finding the k shortest paths, SIAM J. Comput., 1998 - Citeseer.
}


\examples{
data(HIVData)

res<-RchyOptimyx(Signs, -log10(LogRankPvals), '2111012110', 5,FALSE)
plot(res, phenotypeScores=-log10(LogRankPvals), ylab='-log10(Pvalue)')

res<-RchyOptimyx(Signs, OverlapScores, '1201100221', 5,FALSE)
plot(res, phenotypeScores=OverlapScores, ylab='Purity')
}

\author{
  Nima Aghaeepour \email{<naghaeep@bccrc.ca>} and Adrin Jalali \email{<ajalali@bccrc.ca>}
}

\keyword{ FlowCytData }
