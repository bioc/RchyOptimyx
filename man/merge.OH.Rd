\name{merge}
\alias{merge}
\alias{merge-methods}
\alias{merge,OptimizedHierarchy,OptimizedHierarchy-method}
\title{
  Merges two OptimizedHierarchy objects.
}
\description{
  Merges two OptimizedHierarchy objects.
}
\usage{
\S4method{merge}{OptimizedHierarchy,OptimizedHierarchy}(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    First OptimizedHierarchy.
  }
  \item{y}{
    Second OptimizedHierarchy.
  }
}
\value{
  \item{OptH}{An OptimizedHierarchy object}
}


\examples{
data(HIVData)

res1<-RchyOptimyx(Signs, -log10(LogRankPvals), paste(Signs['KI-67+CD4-CCR5+CD127-',], collapse=''), 4,FALSE)
res2<-RchyOptimyx(Signs, -log10(LogRankPvals), paste(Signs['CD45RO-CD8+CD57+CCR5-CD27+CCR7-CD127-',], collapse=''), 2,FALSE)
res3<-RchyOptimyx(Signs, -log10(LogRankPvals), paste(Signs['CD28-CD45RO+CD4-CD57-CD27-CD127-',], collapse=''), 2,FALSE)

res=merge(res1, res2)
res=merge(res,res3)

plot(res, phenotypeScores=-log10(LogRankPvals), ylab='-log10(Pvalue)')
}

\author{
  Nima Aghaeepour \email{<naghaeep@bccrc.ca>} and Adrin Jalali \email{<ajalali@bccrc.ca>}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FlowCytData }
