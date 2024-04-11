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

coverage_barplot <- function(bar_heights, bar_positions, xlim = c(0, 0), ylim = c(0, 0), font_size = 1, width = 0.5){

  if (ylim[2] == 0){
      ylim[2] <- max(bar_heights)
  }
  if (xlim[2] == 0){
      xlim = range(bar_positions)
  }
  
  plot(bar_heights, type="n", xlab="Coverage", ylab="Frequency",
       ylim = ylim, xlim= xlim,
       cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  for ( i in 1:length(bar_heights)){
    rect(bar_positions[i] - width, 0, bar_positions[i] + width, bar_heights[i], col = 'deepskyblue', border = F)
  }
}