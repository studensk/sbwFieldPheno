#' Metadata for development curves
#'
#' @format
#' A data frame with 3000 rows and 4 columns:
#' \describe{
#'   \item{index}{Curve index}
#'   \item{stage}{Development stage; one of L2, L3, L4, L5, L6}
#'   \item{run.index}{Posterior sample index}
#'   \item{colony}{SBW colony membership; one of AB, IPQL, NB, NWT, ON, QC}
#' }
"curve_info"


#' Estimated development rates at treatment temperatures
#'
#' @format
#' A data frame with 27000 rows and 3 columns:
#' \describe{
#'   \item{index}{Curve index}
#'   \item{temp}{Treatment temperature}
#'   \item{rate}{Estimated Development Rate}
#' }
"curves"
