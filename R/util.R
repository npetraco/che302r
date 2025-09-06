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


#' @title       Function to download an example script
#' @description The function will download script that is an example
#'
#' @param script_name Name of the script
#' @param url_head Location of the script (optional). If NULL, default location is used.
#' @param download_location Where to download the script. If NULL, default location is used.
#'
#' @details The function will download an example script. Default url_head is currently:
#'
#' https://raw.githubusercontent.com/npetraco/che302r/refs/heads/main/inst/scripts_bank.
#'
#' Default download_location is currently: tempdir().
#'
#' @return The script should pop up.
#'
#' @examples get_script("XXXXXX.R")
#'
#' @export
get_script <- function(script_name, url_head=NULL, download_location=NULL) {

  if(is.null(url_head)) {
    url_head.loc <- "https://raw.githubusercontent.com/npetraco/che302r/refs/heads/main/inst/scripts_bank"
  } else {
    url_head.loc <- url_head
  }

  if(is.null(download_location)) {
    download_location.loc <- tempdir()
  } else {
    download_location.loc <- download_location
  }

  url_loc <- paste0(url_head.loc, "/", script_name)
  #print(url_loc)

  file_loc <- paste0(download_location.loc,"/",script_name)
  #print(file_loc)

  download.file(url_loc, file_loc)
  print("File downloaded to:")
  print(file_loc)
  file.edit(file_loc)

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
