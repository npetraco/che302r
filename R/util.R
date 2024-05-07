#' Easily update the che302r library by installing the current version from the github site
#'
#' Easily update the che302r library by installing the current version from the github site
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
update_che302r <- function() {
  print("Updating che302r")
  remotes::install_github("npetraco/che302r")
  print("Done!")
}
