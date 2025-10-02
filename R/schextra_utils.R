#' Bivariate binning of single cell data into hexagon cells.
#'
#' \code{.schextra_bin} returns an schex-style binning list.
#'
#' @param obj A SCUBA-supported single-cell object.
#' @param nbins The number of bins partitioning the range of the first
#'    component of the chosen dimension reduction.
#' @param dr A string indicating the reduced dimension
#'    result to calculate hexagon cell representation of.
#' @param use_dims A vector of two integers specifying the dimensions used.
#'
#' @details ...
#'
#' @return A schex-style binning list.
#' @importFrom hexbin hexbin hcell2xy
#' @import SCUBA

.schextra_bin <- function(obj, nbins, dr, use_dims) {
  
    dr <- fetch_reduction(obj, reduction = dr, dims = use_dims)

  # SCUBA only returns the used dims, so this is simplified relative to schex
    xbnds <- range(c(dr[, 1]))
    ybnds <- range(c(dr[, 2]))

    drhex <- hexbin(dr[, 1],
        dr[, 2],
        nbins,
        xbnds = xbnds,
        ybnds = ybnds,
        IDs = TRUE
    )
    cID <- drhex@cID
    drhex <- cbind(
        as.numeric(hcell2xy(drhex)$x),
        as.numeric(hcell2xy(drhex)$y),
        as.numeric(drhex@count)
    )

    colnames(drhex) <- c("x", "y", "number_of_cells")

    res <- list(cID = cID, hexbin.matrix = drhex)

    return(res)
}
