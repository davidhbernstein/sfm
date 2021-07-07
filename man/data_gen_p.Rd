\name{data_gen_p}
\alias{data_gen_p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to generate panel sfa data sets}
\description{This function generates stochastic frontier data sets}
\usage{
data_gen_p(t, N, rand, sig_u, sig_v, sig_r, sig_h, cons, tau, mu, beta1, beta2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ time periods}
  \item{N}{Individuals}
  \item{rand}{seed}
  \item{sig_u}{sigma u}
  \item{sig_v}{sigma v}
  \item{sig_r}{sigma r}
  \item{sig_h}{sigma h}
  \item{cons}{constant}
  \item{tau}{tau of Chen et al 2014}
  \item{mu}{mu of Wang Ho 2010}
  \item{beta1}{beta1}
  \item{beta2}{beta2}
}
\details{NA}
\value{returns pdata.frame object containing sfa variables }
\references{
%% ~put references to the literature/web site here ~
}
\author{David Bernstein}
\note{NA}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{data_gen_p, to see all the true data generating processes}
\examples{

##To install in R: 

library(devtools)
install_github("davidhbernstein/sfm")

library(sfm) 

data_trial <- data_gen_p(t=10,N=100, rand = 100, 
                       sig_u = 1,  sig_v = 0.3, 
                       sig_r = .2, sig_h = .4, 
                       cons = 0.5, tau = 0.5,
                       mu= 0.5, beta1 = 0.5,
                       beta2 = 0.5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
