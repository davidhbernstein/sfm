\name{sfm-package}
\alias{sfm-package}
\alias{sfm}
\docType{package}
\title{
\packageTitle{sfm}}
\description{\packageDescription{sfm}}
\details{The DESCRIPTION file:
\packageDESCRIPTION{sfm}
\packageIndices{sfm}}
\references{Chen, Schmidt, and Wang (2014, JoE); Fillipini and Greene (2016, JPA), Wang and Ho (2010, JoE)}
\keyword{package}
\seealso{\code{\url{http://www.davidharrybernstein.com/software}{http://www.davidharrybernstein.com/software}}}
\examples{
To install in R: 

library(devtools)
install_github("davidhbernstein/sfm", build_vignettes=TRUE)  

## Simple application of the generalized true random effects estimator.

library(sfm)

data_trial <- data_gen_p(t=10,N=100,  rand = 100, 
                         sig_u = 1,   sig_v = 0.3, 
                         sig_r = .2,  sig_h = .4, 
                         cons  = 0.5, beta1 = 0.5,
                         beta2 = 0.5)

psfm(formula    = y_gtre ~ x1 + x2,    
     model_name = "GTRE", 
     data       = data_trial,
     individual = "name",
     PSopt      = FALSE)
}               
                   