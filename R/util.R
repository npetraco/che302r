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
  unloadNamespace("che302r")
  print("Unloaded previous che302r")
  remotes::install_github("npetraco/che302r", force = T)
  print("Installed che302r from github")
  library(che302r)
  print("Loaded updated che302r")
  print("Done!")
}
