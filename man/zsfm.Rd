\name{zsfm}
\alias{zsfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zsfm(formula, model_name = c("ZISF", "ZISF_Z"), data, maxit.bobyqa = 10000, maxit.psoptim = 1000, maxit.optim = 1000, REPORT = 1, trace = 2, pgtol = 0, start_val = FALSE, PSopt = FALSE, bob = TRUE, optHessian = TRUE, inefdec = TRUE, upper = NA, Method = "L-BFGS-B", eta = 0.01, alpha = 0.2, logit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{model_name}{
%%     ~~Describe \code{model_name} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{maxit.bobyqa}{
%%     ~~Describe \code{maxit.bobyqa} here~~
}
  \item{maxit.psoptim}{
%%     ~~Describe \code{maxit.psoptim} here~~
}
  \item{maxit.optim}{
%%     ~~Describe \code{maxit.optim} here~~
}
  \item{REPORT}{
%%     ~~Describe \code{REPORT} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{pgtol}{
%%     ~~Describe \code{pgtol} here~~
}
  \item{start_val}{
%%     ~~Describe \code{start_val} here~~
}
  \item{PSopt}{
%%     ~~Describe \code{PSopt} here~~
}
  \item{bob}{
%%     ~~Describe \code{bob} here~~
}
  \item{optHessian}{
%%     ~~Describe \code{optHessian} here~~
}
  \item{inefdec}{
%%     ~~Describe \code{inefdec} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
  \item{Method}{
%%     ~~Describe \code{Method} here~~
}
  \item{logit}{
%%     ~~Describe \code{logit} here~~
}
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
HERE 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
