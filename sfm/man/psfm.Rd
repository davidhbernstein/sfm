\name{psfm}
\alias{psfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{psfm}
\description{Main panel function of package}
\usage{
psfm(formula, model_name = c("TRE_Z", "GTRE_Z", "TRE", "GTRE", "WMLE", "FD", "GTRE_SEQ1", "GTRE_SEQ2"), data, maxit.bobyqa = 100, maxit.psoptim = 10, maxit.optim = 10, REPORT = 1, trace = 3, pgtol = 0, individual, halton_num = NULL, start_val = FALSE, gamma = FALSE, PSopt = FALSE, bob = TRUE, optHessian = TRUE, inefdec = TRUE, Method = "L-BFGS-B")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description for the model to be estimated}
  \item{model_name}{model name for the estimation}
  \item{data}{a pdata.frame}
  \item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
  \item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
  \item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
  \item{REPORT}{reporting parameter}
  \item{trace}{trace}
  \item{pgtol}{pgtol}
  \item{individual}{individual unit in the regression model}
  \item{halton_num}{number of Halton draws to use in SML models}
  \item{start_val}{starting value (optional)}
  \item{gamma}{gamma}
  \item{PSopt}{use psoptim optimization routine (T/F)}
  \item{bob}{use bobyqa optimization routine (T/F)}
  \item{optHessian}{Logical. Should a numerically differentiated Hessian matrix be returned? (for optim routine)}
  \item{inefdec}{Production or cost function}
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
\references{NA}
\author{
David Bernstein
}
\note{NA}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{NA}
\examples{
data_trial <- data_gen_p(t=10,N=100, rand = 100, 
                         sig_u = 1,  sig_v = 0.3, 
                         sig_r = .2, sig_h = .4, 
                         cons = 0.5, beta1 = 0.5,
                         beta2 = 0.5)

## Trial runs:
max_fd          <-   psfm(formula = y_fd ~  x_fd|z_fd,
                          model_name = "FD",  
                          data=data_trial, 
                          individual = "name")


max_gtre         <-   psfm(formula = y_gtre ~ x1+x2,    
                           model_name = "GTRE", 
                           data=data_trial,
                           individual = "name")

max_tre    <-   psfm(formula = y_tre ~  x1+x2,    
                     model_name = "TRE",  
                     data=data_trial, 
                     individual = "name")

max_wmle   <-   psfm(formula = y_wmle~ x1_w + x2_w,
                     model_name = "WMLE",
                     data=data_trial,
                     individual = "name")
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
