#----------------------------------------------------------------
#' @title Update che302r package
#' @description Easily update the che302r library by installing the current version from the github site.
#'
#' @details Easily update the che302r library by installing the current version from the github site
#' @return Nothing
#'
#' @export
#----------------------------------------------------------------
update_che302r <- function() {
  print("Updating che302r")
  unloadNamespace("che302r")
  print("Unloaded previous che302r")
  remotes::install_github("npetraco/che302r", force = T)
  print("Installed che302r from github")
  library(che302r)
  print("Loaded updated che302r")
  print("Done!")
}


#----------------------------------------------------------------
#' @title Log sum exp trick
#' @description Log sum exp trick to prevent over/underflow
#'
#' @param logv vector of log-ed quantities to add up
#'
#' @details Handy for calculating logZ from a vector of log Boltzmann factors.
#'
#' @return The function will XX
#'
#' @export
#----------------------------------------------------------------
logsumexp <- function(logv)
{
  n <- length(logv)
  max.logv <- max(logv)

  answer <-  max.logv + log(cumsum(c(0,exp(logv - max.logv)))[n+1])

  return(answer)

}
