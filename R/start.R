.onAttach <- function(libname, pkgname) {
  packageStartupMessage("The 'davidhbernstein/sfm' Github package has moved to CRAN with the new name 'sfa'.  
  To install 'sfa' paste:
  install.packages('sfa')
  library(sfa) ")
}