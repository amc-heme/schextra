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
#' @param min_cutoff Minimum cutoff value for feature expression display.
#'    Can be a numeric value (e.g., 0.5) or a quantile string (e.g., "q10" for
#'    10th percentile). Values below this threshold will be capped at the cutoff.
#'    If NULL (default), no minimum cutoff is applied.
#' @param max_cutoff Maximum cutoff value for feature expression display.
#'    Can be a numeric value (e.g., 3) or a quantile string (e.g., "q95" for
#'    95th percentile). Values above this threshold will be capped at the cutoff.
#'    If NULL (default), no maximum cutoff is applied. When split_by is used,
#'    cutoffs are calculated globally across all groups.
#' @param title A string containing the title of the plot. If NULL, defaults
#'    to "Feature (action)".
#' @param xlab A string containing the title of the x axis. If NULL (default),
#'    uses the proper dimension name from the reduction (e.g., "UMAP_1", "PC_1").
#' @param ylab A string containing the title of the y axis. If NULL (default),
#'    uses the proper dimension name from the reduction (e.g., "UMAP_2", "PC_2").
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#'
#' @export
#' @import ggplot2
#' @import rlang
#' @import SCUBA
#' @importFrom dplyr as_tibble
#' @importFrom cowplot theme_cowplot
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
#'
#' # Cap maximum expression at 95th percentile
#' plot_schextra_feature(seurat_obj, feature = "CD8A", max_cutoff = "q95")
#'
#' # Cap maximum at fixed value of 3
#' plot_schextra_feature(seurat_obj, feature = "CD8A", max_cutoff = 3)
#'
#' # Apply both min and max cutoffs
#' plot_schextra_feature(seurat_obj, feature = "CD8A", 
#'                       min_cutoff = "q10", max_cutoff = "q90")
#'
#' # Cutoffs work with split_by (applied globally)
#' plot_schextra_feature(seurat_obj, feature = "CD8A", 
#'                       split_by = "cell_type", max_cutoff = "q95")
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
    min_cutoff = NULL,
    max_cutoff = NULL,
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
    
    # Resolve the feature to expression values. Resolution semantics:
    #   - Keyed names (e.g. "rna_ACTG1") are unambiguous and route directly
    #     to the named assay.
    #   - When `assay` is supplied, only that assay is searched: a bare name
    #     is keyed to force assay-only resolution (bypassing SCUBA's
    #     metadata-first behavior for bare names). A name that is not in that
    #     assay errors -- it is NOT resolved from metadata.
    #   - When `assay` is NULL, resolution is metadata-first then assay
    #     (SCUBA's default for bare names), with a metadata fallback for
    #     names that fetch_feature cannot route to an assay at all.
    # NOTE: fetch_feature() returns a data.frame with cells as ROWS and the
    # feature as a single COLUMN. The returned column name may be keyed
    # (e.g. "rna_ACTG1"), so values are extracted from the returned column
    # rather than by the user-supplied `feature` string.

    # Assay key prefixes (named: assay -> key). Empty if unavailable
    # (e.g. non-Seurat objects), in which case behavior degrades to passing
    # the bare name through to fetch_feature.
    assay_keys <- tryCatch(
        SCUBA::all_keys(obj),
        error = function(e) character(0)
    )

    # Is the user-supplied name already keyed for some assay?
    is_keyed <- length(assay_keys) > 0 &&
        any(vapply(
            assay_keys,
            function(k) nzchar(k) && startsWith(feature, k),
            logical(1)
        ))

    expr_data <- NULL

    if (is_keyed) {
        # Keyed name: unambiguous; route directly (independent of `assay`).
        expr_data <- tryCatch({
            fetch_feature(obj, features = feature, assay = NULL, layer = layer)
        }, error = function(feature_error) {
            stop(sprintf(
                "Failed to retrieve feature '%s': %s",
                feature, conditionMessage(feature_error)
            ))
        })
    } else if (!is.null(assay)) {
        # Assay supplied + bare name: search ONLY that assay. Key the name to
        # bypass SCUBA's metadata-first resolution of bare names.
        key <- if (assay %in% names(assay_keys)) assay_keys[[assay]] else ""
        lookup <- if (nzchar(key)) paste0(key, feature) else feature
        expr_data <- tryCatch({
            fetch_feature(obj, features = lookup, assay = assay, layer = layer)
        }, error = function(feature_error) {
            # Strict assay-only: do NOT fall back to metadata.
            stop(sprintf(
                "Feature '%s' not found in assay '%s'.", feature, assay
            ))
        })
        if (is.null(expr_data) || ncol(expr_data) == 0) {
            stop(sprintf(
                "Feature '%s' not found in assay '%s'.", feature, assay
            ))
        }
    } else {
        # No assay + bare name: metadata-first then assay (SCUBA default),
        # with a metadata fallback when fetch_feature cannot route the name.
        # This preserves the BPCells fast path for real assay features.
        expr_data <- tryCatch({
            fetch_feature(obj, features = feature, assay = assay, layer = layer)
        }, error = function(feature_error) {
            metadata_vars <- tryCatch(
                SCUBA::meta_varnames(obj),
                error = function(e) character(0)
            )
            if (feature %in% metadata_vars) {
                SCUBA::fetch_metadata(obj, vars = feature)
            } else {
                stop(sprintf(
                    "Failed to retrieve feature '%s': %s",
                    feature, conditionMessage(feature_error)
                ))
            }
        })
        # fetch_feature may return a zero-column frame rather than erroring;
        # attempt the metadata fallback before giving up.
        if (is.null(expr_data) || ncol(expr_data) == 0) {
            metadata_vars <- tryCatch(
                SCUBA::meta_varnames(obj),
                error = function(e) character(0)
            )
            if (feature %in% metadata_vars) {
                expr_data <- SCUBA::fetch_metadata(obj, vars = feature)
            }
        }
    }

    # A single feature is requested, so expression lives in one column. The
    # returned column name may be keyed (e.g. "rna_ACTG1"); extract by column
    # rather than by the user-supplied `feature`.
    if (is.null(expr_data) || ncol(expr_data) == 0) {
        stop(sprintf("Feature '%s' not found in the specified assay", feature))
    }

    raw_feature_values <- expr_data[[1]]
    if (!is.numeric(raw_feature_values) && !is.integer(raw_feature_values)) {
        stop(sprintf(
            "Feature '%s' resolved to non-numeric values (class: %s). Only numeric features/metadata can be plotted.",
            feature, paste(class(raw_feature_values), collapse = "/")
        ))
    }
    feature_values <- as.numeric(raw_feature_values)
    
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
    
    # Get proper dimension names from SCUBA
    dim_names <- reduction_dimnames(
        object = obj,
        reduction = dimension_reduction,
        dims = use_dims
    )
    
    # Set plot defaults
    if (is.null(title)) {
        title <- paste0(feature, " (", action, ")")
    }
    
    if (is.null(xlab)) {
        xlab <- dim_names[1]
    }
    
    if (is.null(ylab)) {
        ylab <- dim_names[2]
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
        
        # Apply cutoffs globally across all groups
        out_df$feature_value <- .apply_cutoffs(out_df$feature_value, min_cutoff, max_cutoff)
    } else {
        # No splitting workflow: aggregate feature per bin
        aggregated_values <- .calculate_feature_per_bin(feature_values, out$cID, action)
        
        # Create data frame with bin coordinates and feature values
        # The aggregated_values are in the same order as the bins in hexbin_matrix
        out_df <- as_tibble(out[[2]])
        out_df$feature_value <- aggregated_values
        
        # Apply cutoffs
        out_df$feature_value <- .apply_cutoffs(out_df$feature_value, min_cutoff, max_cutoff)
    }
    
    # Create ggplot with Seurat-style theme
    p <- ggplot(out_df, aes(x = !!sym("x"), y = !!sym("y"), fill = !!sym("feature_value"))) +
        geom_hex(stat = "identity") +
        scale_fill_viridis_c() +
        theme_cowplot(font_size = 14) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(face = "bold")
        ) +
        labs(fill = feature, x = xlab, y = ylab) +
        ggtitle(title)
    
    # Add faceting if split_by is used
    if (!is.null(split_by)) {
        p <- p + facet_wrap(vars(!!sym(split_by)), ncol = ncol, scales = scales)
    }
    
    return(p)
}
