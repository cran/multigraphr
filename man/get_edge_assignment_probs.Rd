% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_edge_assignment_probs.R
\name{get_edge_assignment_probs}
\alias{get_edge_assignment_probs}
\title{Edge assignment probabilities under IEA model given fixed degrees}
\usage{
get_edge_assignment_probs(m, deg.seq, model)
}
\arguments{
\item{m}{integer giving number of edges in multigraph.}

\item{deg.seq}{vector of integers with the sum equal to 2\code{m} representing
the degree sequence of the multigraph.}

\item{model}{character string, either \code{'IEAS'} or \code{'ISA'}.}
}
\value{
A numeric vector representing the edge assignment probabilities
to all possible vertex pair sites. The number of vertex pair sites is given by \eqn{n(n+1)/2}.
}
\description{
Calculates the edge assignment probabilities
given specified degree sequence under the two ways in which the RSM
model can be approximated by the IEA model: \cr
 - the IEAS (independent edge assignment of stubs) model, \cr
 - the ISA (independent stub assignment) model.
}
\details{
The IEAS and ISA edge assignment probabilities to
possible vertex pairs are calculated given a fixed degree sequence \code{deq.seq}
under the IEAS model, and \code{deg.seq}/2\code{m} under the ISA model.

Number of possible vertex pair sites (and thus the length of the edge assignment sequence) is given by
\eqn{(n+1)n/2} where \emph{n} is number of vertices.
}
\examples{
# Under the IEAS model with 10 possible vertex pair sites (4 vertices)
get_edge_assignment_probs(m = 8, deg.seq = c(4, 4, 4, 4), model = "IEAS")

# Under the ISA model with 21 possible vertex pair sites (6 vertices)
get_edge_assignment_probs(m = 10, deg.seq = c(8, 4, 2, 2, 2, 2), model = "ISA")
}
\references{
Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
\cr

Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
\cr

#' Shafie, T., Schoch, D. (2021). Multiplexity analysis of networks using multigraph representations.
\emph{Statistical Methods & Applications} 30, 1425–1444.
}
\seealso{
\code{\link{get_degree_seq}}, \code{\link{iea_model}}
}
\author{
Termeh Shafie
}
