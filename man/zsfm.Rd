\name{zsfm}
\alias{zsfm}
zsfm
\title{
Zero-Inflated Stochastic Frontier Model 
}
\description{
Code to use the Zero-Inflated Stochastic Frontier Model 
}
\usage{
zsfm(formula, 
model_name = c("ZISF", "ZISF_Z"), 
data, 
maxit.bobyqa = 10000,
maxit.psoptim = 1000, 
maxit.optim = 1000, 
REPORT = 1, 
trace = 2, 
pgtol = 0,
start_val = FALSE,
PSopt = FALSE, 
bob = TRUE, 
optHessian = TRUE, 
inefdec = TRUE, 
upper = NA, 
Method = "L-BFGS-B",
eta = 0.01, 
alpha = 0.2, 
logit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{a symbolic description for the model to be estimated}
\item{model_name}{model name for the estimation}
\item{data}{A data set}
\item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
\item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
\item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
\item{REPORT}{reporting parameter}
\item{trace}{trace}
\item{pgtol}{pgtol}
\item{start_val}{starting value (optional)}
\item{PSopt}{use psoptim optimization routine (T/F)}
\item{bob}{use bobyqa optimization routine (T/F)}
\item{optHessian}{Logical. Should a numerically differentiated Hessian matrix be returned? (for optim routine)}
\item{inefdec}{Production or cost function}
\item{upper}{Vector of upper values for the optim package.}
\item{Method}{The method to be used for optim. See ‘Details’ within optim.}
\item{logit}{Choice of using logit function}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chris Parmeter and David Bernstein
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(sfm)

eqz         <- y ~ q1 + q2 + q3 + q4 + q5 + w1 + w2 + w3 + w4|z

data(panel89)

output   <- zsfm(formula    = eqz,
                 model_name = "ZISF_Z",
                 data       = panel89,
                 logit      = TRUE)
}
 