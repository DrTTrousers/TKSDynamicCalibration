# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Basic import namespace populator
#'
#' @return
#' @export
#' @import dplyr
#' @import readr
#' @import DT
#' @import ggplot2
#' @import reshape2
#' @rawNamespace import(plotly, except = last_plot)
#'
#' @examples
#' setup()
setup <- function(){
  print("Imported Dependencies")
}


