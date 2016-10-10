\name{getProjectsList}
\alias{getProjectsList}
\title{Download the EMG projects list into a data frame}

\description{
This downloads an EMG projects list (detailed info) into an R data frame for subsequent analysis.
}

\usage{
getProjectsList()
}

\value{
An R data frame containing a list of all public EMG projects, one row per project.
}

\seealso{\code{\link{getProjectSummary}}, \code{\link{projectSamples}}}

\examples{
pl = getProjectsList()
str(pl)
biggies = pl$Study.ID[pl$Number.Of.Samples >= 50]
ps = getProjectSummary(biggies[3])
}

\keyword{EMG}
