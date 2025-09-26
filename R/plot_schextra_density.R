#' Plot of density of observations from single cell data
#'    in bivariate hexagon cells.
#'
#' @param obj A SCUBA-supported single-cell object.
#' @param nbins The number of bins partitioning the range of the first
#'    component of the chosen dimension reduction.
#' @param dimension_reduction A string indicating the reduced dimension
#'    result to calculate hexagon cell representation of.
#' @param use_dims A vector of two integers specifying the dimensions used.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @import rlang
#' @export
#' 
plot_schextra_density <- function(
    obj,
    nbins = 80,
    dimension_reduction = "UMAP",
    use_dims = c(1,2),
    title = NULL,
    xlab = NULL,
    ylab = NULL) {
  
    out <- .schextra_bin(obj, nbins, dimension_reduction, use_dims)

    if (is.null(title)) {
        title <- "Density"
    }

    if (is.null(xlab)) {
        xlab <- "x"
    }

    if (is.null(ylab)) {
        ylab <- "y"
    }

    out <- as_tibble(out[[2]])

    ggplot(out, aes(x = !!sym("x"), y = !!sym("y"), fill = !!sym("number_of_cells"))) +
        geom_hex(stat = "identity") +
        scale_fill_viridis_c() +
        theme_classic() +
        theme(legend.position = "bottom") +
        ggtitle(title) +
        labs(x = xlab, y = ylab) +
        theme(legend.title = element_blank())
}