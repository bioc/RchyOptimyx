\name{HIVData}
\alias{HIVData}
\alias{OverlapScores}
\alias{LogRankPvals}
\alias{Signs}
\docType{data}
\title{HIVData}
\description{
  A dataset of two sets of scores (particularly, correlation with protection against HIV and overlap with the Naive T-cell population) assigned to immunophenotypes measured by flow cytometry. 10 markers were measured: KI-67, CD28, CD45RO, CD8, CD4, CD57, CCR5, CD27, CCR7, and CD127.
}
\usage{data(HIVData)}
\details{
  This dataset consists of a matrix (\code{Signs}) and 2 numeric vectors
  (\code{LogRankPvals} and \code{OverlapScores}). The \code{Signs}
  matrix consists of 10 columns (one per measured marker) and
  $3^10-1=59048$ rows (one per immunophenotype). See Aghaeepour
  et.al. or the flowType package for more details. For every phenotype (row), the entity
  corresponding to a given marker (column) can be 0, 1, and 2 for
  negative, neutral, and positive respectively. The rownames and
  column names are set respectively.

  \code{LogRankPvals} is a vector of log-rank test P-values to
  determine the correlation between HIV's progression and each of the
  measured immunophenotypes in 466 HIV positive patients (lower values
  represent a stronger correlation). For more details see the Aghaeepour
  et.al. manuscript below.
  The names of the vector match the names of the \code{Signs} matrix. 

  Ganesan et. al. define Naive T-cells as CD28+CD45RO-CD57-CCR5-CD27+
  CCR7+ within the CD3+CD14- compartment. The \code{OverlapScores}
  vector has the proportion of Naive T-cells (as defined above) to the
  total number of cells in any given immunophenotype (a higher value
  represents a larger overlap). The names of the vector match the names of the \code{Signs} matrix. 
}
%\source{
%}

\author{
  Nima Aghaeepour \email{<naghaeep@bccrc.ca>} and Adrin Jalali \email{<ajalali@bccrc.ca>}
}

\references{
  Nima Aghaeepour, Pratip K. Chattopadhyay , Anuradha Ganesan, Kieran O'Neill, Habil Zare, Adrin Jalali, Holger H. Hoos, Mario Roederer, and Ryan R. Brinkman. Early Immunologic Correlates of HIV Protection can be Identified from Computational Analysis of Complex Multivariate T-cell Flow Cytometry Assays. Bioinformatics, pending on minor revisions, 2012.

  Anuradha Ganesan, Pratip K Chattopadhyay, Tess M. Brodie, Jing Qin, Wenjuan Gu, John R. Mascola, Nelson L. Michael, Dean A. Follmann, and Mario Roederer. Immunologic and virologic events in early HIV infection predict subsequent rate of progression. Journal of Infectious Diseases, 2010.
}

