\name{mergeOtu}
\alias{mergeOtu}
\title{Merges two OTU count data data frames to produce a combined OTU count data data frame}

\description{
This function merges two OTU count data data frames, such as returned by \code{getRunOtu}. The merged OTU data is returned as a data frame. The ordering of the two data frames to be merged should not be significant.
}

\usage{
mergeOtu(otu1,otu2)
}

\arguments{
  \item{otu1}{The first OTU data frame to be merged.}
  \item{otu2}{The second OTU data frame to be merged.}
}

\value{
An R data frame containing OTU count data. The first column is the OTU ID, the second is the count, and the third contains a taxonomic classification.
}

\seealso{\code{\link{getSampleOtu}}, \code{\link{getRunOtu}}, \code{\link{projectRuns}}}

\examples{
ps=getProjectSummary("SRP047083")
samp=projectSamples(ps)
runs=runsBySample(ps,samp[2])
otu1=getRunOtu(ps,runs[1])
otu2=getRunOtu(ps,runs[2])
mergeOtu(otu1,otu2)
}

\keyword{EMG}
