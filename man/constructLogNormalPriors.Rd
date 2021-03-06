% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructLogNormalPriors.R
\name{constructLogNormalPriors}
\alias{constructLogNormalPriors}
\title{A function to construct a list of independent priors on covariance parameters from a glmssn object}
\usage{
constructLogNormalPriors(glmssn, std)
}
\arguments{
\item{glmssn}{A fitted glmssn object, preferably with \code{optimOutput$hessian}.}

\item{std}{A numeric vector. This vector specifies the standard deviation of the log-normal priors. See Details for more information.}
}
\value{
A list whose elements are functions parameterised in terms of x, the number of Monte Carlo draws to be taken from the priors when evaluating the expected utility.
}
\description{
The main design function in \code{SSNdesign} is called \code{optimiseSSNDesign}, and it has an argument specifying a list of independent priors on the covariance parameters of a fitted glmssn object.
}
\details{
The argument std can be of length 1, in which case all covariance parameters (parsill, range) will have the same value for standard deviation. Alternatively, it can be of length two, in which case the first element will be used for the partial sill parameters and the second for range parameters. Alternatively, it can be of the same length as glmssn$estimates$theta in which the elements will be matched to their corresponding parameters.
}
\examples{

# Set seed
set.seed(1)

# Create SSN
s <- createSSN(10, binomialDesign(10), 
path = paste(tempdir(), "example00.ssn", sep = "/"), importToR = TRUE)
# And distance matrix
createDistMat(s)

# Simulate data
s <- SimulateOnSSN(s, getSSNdata.frame(s), formula = ~ 1, coefficients = c(1),
 CorModels = c("Spherical.tailup"), CorParms = c(1, 2, 0.1), 
 addfunccol = "addfunccol")$ssn.object

# Fit model
m <- glmssn(s, formula = Sim_Values ~ 1, CorModels = c("Spherical.tailup"),
 addfunccol = "addfunccol")

# Construct log-normal priors
p <- constructLogNormalPriors(m)

}
