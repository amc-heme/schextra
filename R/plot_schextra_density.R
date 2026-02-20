#' Plot of density of observations from single cell data
#'    in bivariate hexagon cells.
#'
#' @param obj A SCUBA-supported single-cell object.
#' @param nbins The number of bins partitioning the range of the first
#'    component of the chosen dimension reduction.
#' @param dimension_reduction A string indicating the reduced dimension
#'    result to calculate hexagon cell representation of.
#' @param use_dims A vector of two integers specifying the dimensions used.
#' @param split_by A string naming a metadata variable to split the plot by.
#'    If NULL (default), no splitting is performed. When provided, creates
#'    separate faceted plots for each unique value in the metadata variable.
#'    Each facet displays the same set of hexagonal bins with cell counts
#'    calculated separately for that metadata group. Bins with zero cells
#'    for a particular group are excluded from that group's facet.
#'    NA values in the metadata will be shown in a separate facet.
#' @param scale_density Logical. If TRUE and split_by is used, each facet
#'    will have an independent color scale showing relative density within that
#'    group (rescaled 0-1). If FALSE (default), all facets share a single color
#'    scale showing absolute cell counts. This allows better visualization of
#'    density patterns in groups with different cell numbers. Only applies when
#'    split_by is used; ignored otherwise.
#' @param ncol Number of columns for facet layout when split_by is used.
#'    If NULL (default), ggplot2 determines the layout automatically.
#' @param scales Should scales be "fixed" (default, same for all facets),
#'    "free" (vary across facets), "free_x", or "free_y"? Only applies
#'    when split_by is used.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis. If NULL (default),
#'    uses the proper dimension name from the reduction (e.g., "UMAP_1", "PC_1").
#' @param ylab A string containing the title of the y axis. If NULL (default),
#'    uses the proper dimension name from the reduction (e.g., "UMAP_2", "PC_2").
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import ggplot2
#' @importFrom dplyr as_tibble group_by mutate ungroup if_else select %>%
#' @import rlang
#' @import SCUBA
#' @importFrom cowplot theme_cowplot
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic density plot
#' plot_schextra_density(seurat_obj, nbins = 80)
#'
#' # Split by cell type with automatic layout
#' plot_schextra_density(seurat_obj, nbins = 80, split_by = "cell_type")
#'
#' # Control layout with 2 columns
#' plot_schextra_density(seurat_obj, split_by = "treatment", ncol = 2)
#'
#' # Free y-axis scales for better contrast per group
#' plot_schextra_density(seurat_obj, split_by = "sample", scales = "free_y")
#'
#' # Use independent color scales per facet for better contrast
#' plot_schextra_density(seurat_obj, split_by = "cell_type", scale_density = TRUE)
#'
#' # Combine independent density scaling with free position scales
#' plot_schextra_density(seurat_obj, split_by = "treatment", 
#'                       scale_density = TRUE, scales = "free_y", ncol = 2)
#' }
#' 
plot_schextra_density <- function(
    obj,
    nbins = 80,
    dimension_reduction = "UMAP",
    use_dims = c(1, 2),
    split_by = NULL,
    scale_density = FALSE,
    ncol = NULL,
    scales = "fixed",
    title = NULL,
    xlab = NULL,
    ylab = NULL) {
  
    # Validate split_by parameter
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
    
    # Validate scale_density parameter
    if (!is.logical(scale_density) || length(scale_density) != 1) {
        stop("scale_density must be a single logical value (TRUE or FALSE)")
    }
    
    # Warn if scale_density is used without split_by
    if (scale_density && is.null(split_by)) {
        warning("scale_density = TRUE has no effect when split_by is NULL. Density scaling only applies to faceted plots.")
    }
  
    out <- .schextra_bin(obj, nbins, dimension_reduction, use_dims)

    # Get proper dimension names from SCUBA
    dim_names <- reduction_dimnames(
        object = obj,
        reduction = dimension_reduction,
        dims = use_dims
    )

    if (is.null(title)) {
        title <- "Density"
    }

    if (is.null(xlab)) {
        xlab <- dim_names[1]
    }

    if (is.null(ylab)) {
        ylab <- dim_names[2]
    }

    # Add split metadata if requested
    if (!is.null(split_by)) {
        # Use expanded data structure with per-group counts
        out_df <- .add_metadata_to_bins(
            obj, 
            out$cID,
            out[[2]],
            out$cell,
            split_by, 
            dimension_reduction, 
            use_dims
        )
    } else {
        # Convert to tibble (no splitting)
        out_df <- as_tibble(out[[2]])
    }
    
    # Apply independent density scaling if requested
    if (!is.null(split_by) && scale_density) {
        # Rescale density values to 0-1 within each group
        # This allows each facet to have its own color scale
        out_df <- out_df %>%
            group_by(!!sym(split_by)) %>%
            mutate(
                min_cells = min(.data$number_of_cells),
                max_cells = max(.data$number_of_cells),
                scaled_density = if_else(
                    .data$max_cells == .data$min_cells,
                    0.5,  # If all bins have same count, use middle value
                    (.data$number_of_cells - .data$min_cells) / (.data$max_cells - .data$min_cells)
                )
            ) %>%
            ungroup() %>%
            select(-.data$min_cells, -.data$max_cells)  # Remove temporary columns
        
        # Use scaled_density for plotting
        fill_var <- "scaled_density"
        legend_label <- "Scaled Density"
    } else {
        # Use raw number_of_cells for plotting
        fill_var <- "number_of_cells"
        legend_label <- NULL  # Will use default (no title with element_blank())
    }

    # Create base plot with Seurat-style theme
    p <- ggplot(out_df, aes(x = !!sym("x"), y = !!sym("y"), fill = !!sym(fill_var))) +
        geom_hex(stat = "identity") +
        scale_fill_viridis_c() +
        theme_cowplot(font_size = 14) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(face = "bold")
        ) +
        ggtitle(title) +
        labs(x = xlab, y = ylab)
    
    # Set legend title based on whether density scaling is used
    if (!is.null(legend_label)) {
        p <- p + labs(fill = legend_label)
    } else {
        p <- p + theme(legend.title = element_blank())
    }

    # Add faceting if split_by is used
    if (!is.null(split_by)) {
        p <- p + facet_wrap(vars(!!sym(split_by)), ncol = ncol, scales = scales)
    }

    return(p)
}