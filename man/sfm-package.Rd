\name{sfm-package}
\alias{sfm-package}
\alias{sfm}
\docType{package}
\title{
\packageTitle{sfm}
}
\description{
\packageDescription{sfm}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{sfm}
\packageIndices{sfm}
This package is currently being developed to work in a wide range of applications.   The main function of the package is "sfm" which can apply four different panel based SFA estimators.
}
\author{
\packageAuthor{sfm}

Maintainer: \packageMaintainer{sfm}
}
\references{Chen Schmidt and Wang 2014 \\ Fillipini and Greene 2016 \\ Wang Ho 2010
}
\keyword{ package }
\seealso{\code{\url{http://www.davidharrybernstein.com/software}{http://www.davidharrybernstein.com/software}}}
\examples{
## To INSTALL: 1) setwd to location of downloaded file -- 2)ensure not extra folder
## 3) install devtools and roxygen2, 4) install other packages bellow and 5) install sfm: 

## be sure all these packages are installed for sfm to run  
list.of.packages = c("plm","minqa","randtoolbox","matrixStats","frontier","Jmisc","mnormt", "devtools", "roxygen2") 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {install.packages(new.packages)}
lapply(list.of.packages, require, character.only=T)  
rm(list.of.packages, new.packages)

install("sfm")



## Simple application of the generalized true random effects estimator.

library(sfm)

data("data_trial")

max_gtre         <-   psfm(formula = y_gtre ~ x1+x2,    
                            model = "GTRE", 
                            data=data_trial,
                            individual = "name")