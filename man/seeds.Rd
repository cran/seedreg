% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seeds_function.R
\name{seeds}
\alias{seeds}
\title{Param: Seeds}
\usage{
seeds(data, trat, nrep, time)
}
\arguments{
\item{data}{Data.frame containing the responses of the evaluations in separate columns side by side and without the columns with the identification of the factors}

\item{trat}{Vector of treatments with n repetitions}

\item{nrep}{Number of repetitions}

\item{time}{Vector containing time}
}
\value{
Returns a data.frame with the indices
}
\description{
Simplification of functions: acc, iv, tm and tml.
}
\examples{
data("substrate")
seeds(substrate[,c(3:18)],
      trat = substrate$Trat,
      nrep = 4,
      time = 1:16)
}
