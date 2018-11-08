runShiny <- function() {
  appDir <- system.file("shinyapp", package = "simDeNet")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `simDeNet`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}