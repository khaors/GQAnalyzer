#' @title
#' GQAnalyzer_gui
#' @description
#' GQAnalyzer GUI
#' @import shiny
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
GQAnalyzer_gui <- function() {
  appDir <- system.file("Shiny", "GQAnalyzer_gui", package = "GQAnalyzer")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `GQAnalyzer`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
