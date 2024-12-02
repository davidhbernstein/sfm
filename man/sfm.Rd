\name{sfm}
\alias{sfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sfm
}
\description{
Cross-sectional stochastic frontier models are estimated with the sfm() call.  For panel models, see the psfm() call.   
}
\usage{
sfm(formula, model_name = c("NHN", "NHN_Z", "NE", "NE_Z", "THT", "NTN", "NR"), data, maxit.bobyqa = 10000, maxit.psoptim = 1000, maxit.optim = 1000, REPORT = 1, trace = 2, pgtol = 0, start_val = FALSE, PSopt = FALSE, bob = TRUE, optHessian = TRUE, inefdec = TRUE, upper = NA, Method = "L-BFGS-B")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{a symbolic description for the model to be estimated}
\item{model_name}{model name for the estimation includes the: normal-half normal (NHN), normal-exponential (NE), student's t-half t (THT),  Normal-Rayleigh (NR), and the normal-truncated normal (NTN).}
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
NA
}
\author{
David Bernstein
}
\note{The options include the normal-half normal (NHN), Normal-exponential (NE), student's t-half t (THT), and the normal-truncated normal (NTN).  NHN_Z and NE_Z are extensions for the NHN and NE models that allow for modeling the u-component of those models with determinants of inefficiency.

Outputs include E[exp(-u)|e] given by exp_u_hat, following Battese and Coelli (1988, JoE), where appropriate.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
cs_data_trial   <- data_gen_cs(N= 1000, rand   = 1,  sig_u  = 0.3, sig_v  = 0.3, cons   = 0.5,       beta1  = 0.5,   beta2  = 0.5, a      = 4,
mu     = 1)

cs.nhnz     <-  sfm(formula    = y_pcs_z ~ x1 +x2| z,    model_name = "NHN",                  
                    data       = cs_data_trial,          PSopt      = TRUE)

cs.nez      <-  sfm(formula    = y_pcs_z ~ x1 +x2| z,    model_name = "NE",                 
                    data       = cs_data_trial,          PSopt      = TRUE)

cs.nhn       <- sfm(formula    = y_pcs ~ x1 +x2,         model_name = "NHN", 
                    data       = cs_data_trial,          PSopt      = TRUE)

cs.ne        <- sfm(formula    = y_pcs_e ~ x1 +x2,       model_name = "NE", 
                    data       = cs_data_trial)

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
