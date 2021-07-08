\name{psfm}
\alias{psfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Main panel function of package}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
psfm(formula, model_name = c("GTRE", "TRE", "PCS", "WMLE", "FD", "GTRE_SEQ1", "GTRE_SEQ2"), data, maxit = 100, REPORT = 1, trace = 3, pgtol = 0, individual, start_val = FALSE, inefdec = TRUE, gamma = FALSE, z_vars, bob = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description for the model to be estimated}
  \item{model_name}{model name}
  \item{data}{a pdata.frame}
  \item{maxit}{Maximum number of iterations}
  \item{REPORT}{reporting}
  \item{trace}{trace}
  \item{pgtol}{pgtol}
  \item{individual}{individual}
  \item{start_val}{starting value}
  \item{inefdec}{Production or cost function}
  \item{gamma}{gamma}
  \item{bob}{optimization routine}
}
\references{
HI
}
\author{
David Bernstein
}
\note{
HI
}

\examples{## R Package
##To install in R: 

library(devtools)
install_github("davidhbernstein/sfm")  

rm(list = ls())

library(sfm)

data_trial <- data_gen_p(t=10,N=100, rand = 100, 
                         sig_u = 1,  sig_v = 0.3, 
                         sig_r = .2, sig_h = .4, 
                         cons = 0.5, tau = 0.5,
                         mu= 0.5, beta1 = 0.5,
                         beta2 = 0.5)

## Trial runs:
max_fd          <-   psfm(formula = y_fd ~  x_fd,
                          z_vars  = c("z_fd"),
                          model_name = "FD",  
                          data=data_trial, 
                          individual = "name")


max_gtre         <-   psfm(formula = y_gtre ~ x1+x2,    
                            model_name = "GTRE", 
                            data=data_trial,
                            individual = "name")

print(cor( max_gtre[[7]] ,exp(-data_trial$u)))


max_tre    <-   psfm(formula = y_tre ~  x1+x2,    
                            model_name = "TRE",  
                            data=data_trial, 
                            individual = "name")

print(cor(max_tre[[7]],exp(-data_trial$u)))

max_pcs    <-   psfm(formula = y_pcs~ x1+x2,
                            model_name = "PCS",
                            data=data_trial,
                            individual = "name")

print(cor(max_pcs[[7]],exp(-data_trial$u) ))


max_wmle   <-   psfm(formula = y_wmle~ x1_w + x2_w,
                            model_name = "WMLE",
                            data=data_trial,
                            individual = "name")

print(cor(max_wmle[[6]],exp(-data_trial$u)))

effs <- cbind(max_pcs[[7]], max_tre[[7]],max_wmle[[6]], max_gtre[[7]], exp(-data_trial$u) )

cor_mat_ken  <- cor(effs, method = "kendall")
cor_mat_pea  <- cor(effs, method = "pearson")
cor_mat_spe  <- cor(effs, method = "spearman")

rownames(cor_mat_ken) <- c("PCS", "TRE", "WMLE","GTRE", "True")
rownames(cor_mat_pea) <- c("PCS", "TRE", "WMLE","GTRE", "True")
rownames(cor_mat_spe) <- c("PCS", "TRE", "WMLE","GTRE", "True")
colnames(cor_mat_ken) <- c("PCS", "TRE", "WMLE","GTRE", "True")
colnames(cor_mat_pea) <- c("PCS", "TRE", "WMLE","GTRE", "True")
colnames(cor_mat_spe) <- c("PCS", "TRE", "WMLE","GTRE", "True")
colnames(effs)        <- c("PCS", "TRE", "WMLE","GTRE", "True")

cor_mat_ken
cor_mat_pea
cor_mat_spe
summary(effs)
}
