\name{Forest}
\alias{Forest}
\docType{data}
\title{
Forest
}
\description{
Cross-sectional data for 3,249 active Norwegian forest owners in 2003, compiled by Statistics Norway from survey and administrative sources. The data include timber output, labour, harvested area, capital input, and property characteristics, and are used in Lien, Størdal and Baardsen (2007) to study technical efficiency in timber production.
\usage{data("Forest")}
\format{
  A data frame with 3249 observations on the following 4 variables.
  \describe{
    \item{\code{Y}}{a numeric vector containing timber harvesting level (\eqn{m^3})}
    \item{\code{X1}}{a numeric vector containing labour input (hours)}
    \item{\code{X2}}{a numeric vector containing forest area cut (ha)}
    \item{\code{X3}}{a numeric vector containing the capital value of maximum sustainable yield (Norwegian kroner)}
  }
}
\details{
Output is measured as annual timber sales, with timber inventories ignored by assuming opening and closing inventories are equal. Labour input includes owner, family, hired, and contractor labour; contractor hours are imputed from reported contractor costs using a fixed hourly rate of 1,925 NOK. Harvested forest area measures hectares subject to final fellings, while capital input is measured as the value of forest increment, calculated as the weighted mean timber price multiplied by the maximum sustainable yield of the property. For a small number of properties with missing timber price or quantity information, values are replaced with sample means, and the sample over-represents relatively large forest properties.
}
\source{
Lien, G., Størdal, S. and Baardsen, S. (2007) ‘Technical efficiency in timber production and effects of other income sources’, Small-scale Forestry, 6(1), pp. 65–78. doi:10.1007/s11842-007-9005-3.
}
\references{
Statistics Norway (2006) ‘Sample Survey of Forestry 2004’, Official Statistics of Norway, D 341. Oslo: Statistics Norway. Available at: https://www.ssb.no/a/publikasjoner/pdf/nos_d341/nos_d341.pdf
}
\examples{
data(Forest)
plot(Forest)
}
\keyword{datasets}
