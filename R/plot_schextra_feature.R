#' Plot feature expression in hexagonal bins
#'
#' Visualizes gene/feature expression across hexagonal bins of dimension
#' reduction space. Expression values are aggregated within each bin using
#' the specified method (mean, median, or sum). On-the-fly binning is
#' performed using the schextra approach, with expression data accessed
#' via SCUBA's fetch_feature() function.
#'
#' @param obj A SCUBA-supported single-cell object.
#' @param feature A string naming the feature/gene to plot.
#' @param assay Name of assay to use. If NULL (default), uses the object's
#'    default/active assay.
#' @param layer Name of assay layer to use (e.g., "data", "counts", "scale.data").
#'    If NULL (default), uses the object's default layer for the specified assay.
#'    For Seurat objects, common options include "data" (normalized), "counts" 
#'    (raw counts), or "scale.data" (scaled).
#' @param action Aggregation method for expression within bins. One of:
#'    "mean" (default), "median", or "sum".
#' @param nbins The number of bins partitioning the range of the first
#'    component of the chosen dimension reduction.
#' @param dimension_reduction A string indicating the reduced dimension
#'    result to use (e.g., "UMAP", "PCA").
#' @param use_dims A vector of two integers specifying the dimensions to plot.
#' @param split_by A string naming a metadata variable to split the plot by.
#'    If NULL (default), no splitting is performed. When provided, creates
#'    separate faceted plots for each unique value in the metadata variable.
#'    Expression is calculated separately for each group within each bin.
#'    Bins with zero cells for a particular group are excluded from that
#'    group's facet. NA values in the metadata will be shown in a separate facet.
#' @param ncol Number of columns for facet layout when split_by is used.
#'    If NULL (default), ggplot2 determines the layout automatically.
#' @param scales Should scales be "fixed" (default, same for all facets),
#'    "free" (vary across facets), "free_x", or "free_y"? Only applies
#'    when split_by is used.
#' @param title A string containing the title of the plot. If NULL, defaults
#'    to "Feature (action)".
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#'
#' @export
#' @import ggplot2
#' @import rlang
#' @import SCUBA
#' @importFrom dplyr as_tibble
#'
#' @examples
#' \dontrun{
#' # Basic feature plot with mean aggregation
#' plot_schextra_feature(seurat_obj, feature = "CD8A")
#'
#' # Use median aggregation
#' plot_schextra_feature(seurat_obj, feature = "CD8A", action = "median")
#'
#' # Plot on PCA instead of UMAP
#' plot_schextra_feature(seurat_obj, feature = "CD8A", 
#'                       dimension_reduction = "PCA")
#'
#' # Use raw counts instead of normalized data
#' plot_schextra_feature(seurat_obj, feature = "CD8A", layer = "counts")
#'
#' # Split by cell type
#' plot_schextra_feature(seurat_obj, feature = "CD8A", split_by = "cell_type")
#'
#' # Custom layout and scales
#' plot_schextra_feature(seurat_obj, feature = "CD8A",
#'                       split_by = "treatment", ncol = 2, scales = "free_y")
#' }

plot_schextra_feature <- function(
    obj,
    feature,
    assay = NULL,
    layer = NULL,
    action = "mean",
    nbins = 80,
    dimension_reduction = "UMAP",
    use_dims = c(1, 2),
    split_by = NULL,
    ncol = NULL,
    scales = "fixed",
    title = NULL,
    xlab = NULL,
    ylab = NULL) {
    
    # Validate action parameter
    valid_actions <- c("mean", "median", "sum")
    if (!action %in% valid_actions) {
        stop(sprintf("action must be one of: %s", paste(valid_actions, collapse = ", ")))
    }
    
    # Validate split_by parameter if provided
    if (!is.null(split_by)) {
        # Validate scales parameter
        valid_scales <- c("fixed", "free", "free_x", "free_y")
        if (!scales %in% valid_scales) {
            stop(sprintf("scales must be one of: %s", paste(valid_scales, collapse = ", ")))
        }
        
        # Check if metadata exists
        tryCatch({
            test_meta <- fetch_metadata(obj, vars = split_by)
        }, error = function(e) {
            available_vars <- tryCatch({
                names(fetch_metadata(obj, full_table = TRUE))
            }, error = function(e2) {
                "unable to retrieve"
            })
            stop(sprintf("split_by variable '%s' not found in object metadata. Available variables: %s",
                         split_by, 
                         paste(available_vars, collapse = ", ")))
        })
        
        # Warn if too many groups
        n_groups <- length(unique(test_meta[[split_by]]))
        if (n_groups > 10) {
            warning(sprintf("split_by variable '%s' has %d unique values. Plot may be difficult to read with many facets.",
                           split_by, n_groups))
        }
    }
    
    # IMPORTANT: fetch_feature() returns a data.frame with:
    #   - Cells as ROWS (one row per cell)
    #   - Features as COLUMNS (one column per feature)
    # This is different from traditional expression matrices where features are rows!
    
    # Get feature expression data using SCUBA
    expr_data <- tryCatch({
        fetch_feature(obj, features = feature, assay = assay, layer = layer)
    }, error = function(e) {
        stop(sprintf("Failed to retrieve feature '%s': %s", feature, e$message))
    })
    
    # Extract feature expression as numeric vector (fetch_feature returns data.frame: cells as rows, features as columns)
    if (ncol(expr_data) == 0 || !(feature %in% colnames(expr_data))) {
        stop(sprintf("Feature '%s' not found in the specified assay", feature))
    }
    
    feature_values <- as.numeric(expr_data[[feature]])
    
    # Validate feature_values length before proceeding
    if (length(feature_values) == 0) {
        stop(sprintf("No expression values retrieved for feature '%s'. The feature may not be present in any cells.", feature))
    }
    
    # Create hexagonal bins
    out <- .schextra_bin(obj, nbins, dimension_reduction, use_dims)
    
    # Validate that feature_values length matches number of cells in binning
    if (length(feature_values) != length(out$cID)) {
        stop(sprintf(
            "Mismatch between feature data and dimension reduction: feature data has %d cells but dimension reduction has %d cells. Ensure the same cells are present in both.",
            length(feature_values),
            length(out$cID)
        ))
    }
    
    # Set plot defaults
    if (is.null(title)) {
        title <- paste0(feature, " (", action, ")")
    }
    
    if (is.null(xlab)) {
        xlab <- "x"
    }
    
    if (is.null(ylab)) {
        ylab <- "y"
    }
    
    # Handle splitting logic
    if (!is.null(split_by)) {
        # Split-by workflow: expand bins by metadata groups with feature values
        out_df <- .add_metadata_to_bins_with_feature(
            obj,
            out$cID,
            out[[2]],
            out$cell,
            feature_values,
            action,
            split_by,
            dimension_reduction,
            use_dims
        )
    } else {
        # No splitting workflow: aggregate feature per bin
        aggregated_values <- .calculate_feature_per_bin(feature_values, out$cID, action)
        
        # Create data frame with bin coordinates and feature values
        # The aggregated_values are in the same order as the bins in hexbin_matrix
        out_df <- as_tibble(out[[2]])
        out_df$feature_value <- aggregated_values
    }
    
    # Create ggplot
    p <- ggplot(out_df, aes(x = !!sym("x"), y = !!sym("y"), fill = !!sym("feature_value"))) +
        geom_hex(stat = "identity") +
        scale_fill_viridis_c() +
        theme_classic() +
        theme(legend.position = "bottom") +
        labs(fill = feature, x = xlab, y = ylab) +
        ggtitle(title)
    
    # Add faceting if split_by is used
    if (!is.null(split_by)) {
        p <- p + facet_wrap(vars(!!sym(split_by)), ncol = ncol, scales = scales)
    }
    
    return(p)
}
