#' coverage_barplot - function to plot a barplot in the "deepskyblue" color.
#'
#' @param bar_heights, a vector with heights of individual bars (defining xlim for now)
#'
#' @param bar_positions, locations where bars are located (defining xlim)
#'
#' @param font_size, size of font (default = 1)
#'
#' @param width, width of individual bars (Default well suited for coverage is 0.5)
#'
#' @export

coverage_barplot <- function(bar_heights, bar_positions, font_size = 1, width = 0.5){

  plot(bar_heights, type="n", xlab="Coverage", ylab="Frequency",
       ylim=c(0, max(bar_heights)), xlim=range(bar_positions),
       cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  for ( i in 1:length(bar_heights)){
    rect(bar_positions[i] - width, 0, bar_positions[i] + width, bar_heights[i], col = 'deepskyblue', border = F)
  }
}
