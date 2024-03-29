% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gof_stats.R
\name{gof_stats}
\alias{gof_stats}
\title{Exact probability distributions and moments of goodness of fit statistics}
\usage{
gof_stats(m, dof, m.seq, prob.mg, Q.seq)
}
\arguments{
\item{m}{integer giving number of edges in multigraph.}

\item{dof}{integer giving degrees of freedom of test performed.}

\item{m.seq}{matrix of integers, each row representing the
edge multiplicity sequence of a multigraph (which correspond to observed values).}

\item{prob.mg}{numerical vector representing a given probability distribution of
multigraphs/edge multiplicity sequences in \code{m.seq}.}

\item{Q.seq}{a numeric vector representing the hypothetical edge assignment probabilities
to all possible vertex pair sites (from which expected values are calculate).}
}
\value{
\item{test.summary}{Expected value and variances of test statistics (\code{stat}),
 critical values (\code{cv}) according to asymptotic chi^2-distribution and
 according to cdf's of test statistics,
 significance level (alpha) according to asymptotic chi^2 distribution,
 power of tests (\code{P(stat>cv)}), critical values and power
 according to the distributions of test statistics (\code{cv(stat)}
 and \code{ P(Stat>cv(Stat))}).}
 \item{degrees.of.freedom}{Degrees of freedom for tests performed.}
 \item{probS}{Probability distributions of Pearson statistic \code{S}.}
 \item{probA}{Probability distributions of information divergence statistic \code{A}.}
 \item{adjusted.stats}{Expected values and variances for adjusted test statistics,
 preferred adjusted statistics.}
 \item{adjusted.chi2}{Degrees of freedom for adjusted  chi^2-distribution.}
 \item{power.apx}{Power approximations according to adjusted statistics.}
}
\description{
Goodness of fit between two specified edge multiplicity sequences
(e.g. observed vs. expected).
Pearson (\emph{S}) and information divergence (\emph{A}) tests statistics are used and
the exact distribution of these statistics,  their asymptotic chi^2-distributions,
and their first two central moments are calculated using this function.
Only practical for small multigraphs.
}
\details{
The tests are performed using goodness-of-fit measures between two edge multiplicity sequences
(e.g. observed vs. expected).

Test statistics of Pearson type (\emph{S}) and
of information divergence (\emph{A}) type are used and summary
of tests given these two statistics are given as output. The adjusted statistics and
chi^2-distributions are useful for better power calculations.
}
\examples{
# Generate a set of edge multiplicity sequences (random multigraphs) and
# its probability distribution using rsm_model() with degree sequence [4,4,6,6]
rsm <- rsm_model(deg.seq = c(4,4,6,6))
mg <- as.matrix(rsm$m.seq)
mg.p <- rsm$prob.dists[, 1]

# Generate edge assignment probabilities from which the second set of
# edge multiplicity sequences is generated from using the iea_model()
deg.f <- (4*5)/2 - 1
eap <- get_edge_assignment_probs(m = 10,
                   deg.seq = c(4,4,6,6), model = 'IEAS')

# Perform the test
test <- gof_stats(m = 10, dof = deg.f,
                   m.seq = mg, prob.mg = mg.p, eap)

}
\references{
Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
\cr

Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
}
\seealso{
\code{\link{gof_sim}},\code{\link{get_edge_assignment_probs}},
\code{\link{nsumk}}
}
\author{
Termeh Shafie
}
