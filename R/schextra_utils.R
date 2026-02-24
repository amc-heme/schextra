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
#' @return A schex-style binning list with components:
#'   \describe{
#'     \item{cID}{Vector of bin assignments for each cell}
#'     \item{hexbin.matrix}{Matrix with x, y coordinates and cell counts}
#'     \item{cell}{Vector of unique bin IDs (for mapping to matrix rows)}
#'   }
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
    cell <- drhex@cell  # Store bin IDs for mapping
    drhex <- cbind(
        as.numeric(hcell2xy(drhex)$x),
        as.numeric(hcell2xy(drhex)$y),
        as.numeric(drhex@count)
    )

    colnames(drhex) <- c("x", "y", "number_of_cells")

    res <- list(cID = cID, hexbin.matrix = drhex, cell = cell)

    return(res)
}

#' Expand hexbin results by metadata groups
#'
#' Creates an expanded data frame where each bin appears once per metadata
#' group, with cell counts calculated separately for each group. Bins with
#' zero cells for a particular group are excluded from that group's facet.
#' Cells with NA metadata are treated as a separate "NA" group.
#'
#' @param obj A SCUBA-supported single-cell object
#' @param cID Vector of cell-to-bin assignments from hexbin
#' @param hexbin_matrix Matrix with x, y coordinates and total cell counts
#' @param bin_ids Vector of unique bin IDs from hexbin object (used to map bin IDs to matrix rows)
#' @param metadata_var String naming the metadata variable
#' @param dr_name Name of dimension reduction (for cell ordering)
#' @param use_dims Dimensions used (for cell ordering)
#'
#' @return Data frame with columns: x, y, number_of_cells, and metadata_var
#' @import SCUBA
#' @importFrom stats setNames
#' @keywords internal

.add_metadata_to_bins <- function(obj, cID, hexbin_matrix, bin_ids, metadata_var, dr_name, use_dims) {
    
    # Get cell IDs in the correct order
    dr_data <- fetch_reduction(obj, reduction = dr_name, dims = use_dims)
    cell_ids <- rownames(dr_data)
    
    # Fetch metadata for this variable
    metadata_df <- fetch_metadata(obj, vars = metadata_var)
    
    # Create cell_id -> metadata mapping
    cell_metadata <- setNames(metadata_df[[metadata_var]], rownames(metadata_df))
    
    # Match metadata to cells in the order they appear in the reduction
    cell_metadata_ordered <- cell_metadata[cell_ids]
    
    # Store original factor levels if this is a factor (for restoration later)
    is_factor <- is.factor(cell_metadata_ordered)
    if (is_factor) {
        original_levels <- levels(cell_metadata_ordered)
    }
    
    # Convert NAs to "NA" string for proper faceting (only if NAs actually exist)
    has_na <- any(is.na(cell_metadata_ordered))
    if (has_na) {
        if (is_factor) {
            levels(cell_metadata_ordered) <- c(levels(cell_metadata_ordered), "NA")
            original_levels <- levels(cell_metadata_ordered)  # Update stored levels
        }
        cell_metadata_ordered[is.na(cell_metadata_ordered)] <- "NA"
    }
    
    # Get unique bins that have at least one cell
    unique_bins <- unique(cID)
    
    # Create mapping from bin_id to row_index in hexbin_matrix
    # bin_ids contains the bin IDs in the same order as hexbin_matrix rows
    bin_id_to_row <- setNames(seq_along(bin_ids), as.character(bin_ids))
    
    # Create expanded data structure
    expanded_rows <- list()
    row_index <- 1
    
    for (bin_num in unique_bins) {
        # Find cells in this bin
        cells_in_bin <- which(cID == bin_num)
        
        # Get metadata values for cells in this bin
        meta_values <- cell_metadata_ordered[cells_in_bin]
        
        # Get unique groups in this bin (only groups actually present)
        # Note: Using unique() instead of table() to avoid including unused factor levels
        unique_groups <- unique(meta_values)
        
        # Get row index for this bin ID
        row_idx <- bin_id_to_row[as.character(bin_num)]
        
        # Get bin coordinates using the row index
        bin_x <- hexbin_matrix[row_idx, "x"]
        bin_y <- hexbin_matrix[row_idx, "y"]
        
        # Create one row per group that has cells in this bin
        for (group in unique_groups) {
            # Count cells for this specific group in this bin
            group_count <- sum(meta_values == group)
            
            expanded_rows[[row_index]] <- data.frame(
                x = bin_x,
                y = bin_y,
                number_of_cells = group_count,
                metadata_group = as.character(group),  # Convert to character to avoid factor level issues
                stringsAsFactors = FALSE
            )
            row_index <- row_index + 1
        }
    }
    
    # Combine all rows into a single data frame
    result <- do.call(rbind, expanded_rows)
    
    # Rename the metadata column to the actual variable name
    colnames(result)[colnames(result) == "metadata_group"] <- metadata_var
    
    # Restore factor type with original levels if input was a factor
    if (is_factor) {
        result[[metadata_var]] <- factor(result[[metadata_var]], levels = original_levels)
    }
    
    return(result)
}

#' Apply minimum and maximum cutoffs to feature values
#'
#' Caps feature values at specified minimum and maximum thresholds.
#' Supports both numeric cutoffs and quantile-based cutoffs (e.g., "q95").
#' When both cutoffs are provided, minimum cutoff is applied first, then
#' maximum cutoff.
#'
#' @param values Numeric vector of feature values
#' @param min_cutoff Minimum cutoff: NULL (no cutoff), numeric value, or
#'    quantile string "q##" where ## is 0-100 (e.g., "q5" for 5th percentile)
#' @param max_cutoff Maximum cutoff: NULL (no cutoff), numeric value, or
#'    quantile string "q##" where ## is 0-100 (e.g., "q95" for 95th percentile)
#'
#' @return Numeric vector with cutoffs applied
#' @importFrom stats quantile
#' @keywords internal

.apply_cutoffs <- function(values, min_cutoff = NULL, max_cutoff = NULL) {
    result <- values
    
    # Apply minimum cutoff
    if (!is.null(min_cutoff)) {
        if (is.character(min_cutoff) && grepl("^q[0-9]+$", min_cutoff)) {
            percentile <- as.numeric(sub("^q", "", min_cutoff))
            if (percentile < 0 || percentile > 100) {
                stop("min_cutoff quantile must be between 0 and 100 (e.g., 'q5')")
            }
            cutoff_value <- quantile(values, probs = percentile / 100, na.rm = TRUE)
        } else if (is.numeric(min_cutoff)) {
            cutoff_value <- min_cutoff
        } else {
            stop("min_cutoff must be NULL, a numeric value, or a quantile string (e.g., 'q5')")
        }
        result <- pmax(result, cutoff_value, na.rm = TRUE)
    }
    
    # Apply maximum cutoff
    if (!is.null(max_cutoff)) {
        if (is.character(max_cutoff) && grepl("^q[0-9]+$", max_cutoff)) {
            percentile <- as.numeric(sub("^q", "", max_cutoff))
            if (percentile < 0 || percentile > 100) {
                stop("max_cutoff quantile must be between 0 and 100 (e.g., 'q95')")
            }
            cutoff_value <- quantile(values, probs = percentile / 100, na.rm = TRUE)
        } else if (is.numeric(max_cutoff)) {
            cutoff_value <- max_cutoff
        } else {
            stop("max_cutoff must be NULL, a numeric value, or a quantile string (e.g., 'q95')")
        }
        result <- pmin(result, cutoff_value, na.rm = TRUE)
    }
    
    return(result)
}

#' Aggregate feature expression values within hexagonal bins
#'
#' Calculates aggregated feature expression for each hexagonal bin using
#' the specified aggregation method (mean, median, or sum).
#'
#' @param feature_values Numeric vector of expression values (one per cell)
#' @param cID Vector of cell-to-bin assignments from hexbin
#' @param action Aggregation method: "mean", "median", or "sum"
#'
#' @return Named numeric vector where names are bin IDs and values are
#'    aggregated expression levels for each bin
#' @importFrom stats median
#' @keywords internal

.calculate_feature_per_bin <- function(feature_values, cID, action) {
    
    # Validate inputs
    if (!is.numeric(feature_values)) {
        stop("feature_values must be numeric")
    }
    
    if (!action %in% c("mean", "median", "sum")) {
        stop("action must be one of: 'mean', 'median', 'sum'")
    }
    
    # Define aggregation function based on action
    agg_func <- switch(action,
        mean = function(x) mean(x, na.rm = TRUE),
        median = function(x) median(x, na.rm = TRUE),
        sum = function(x) sum(x, na.rm = TRUE)
    )
    
    # Aggregate expression values by bin
    aggregated <- tapply(feature_values, cID, FUN = agg_func)
    
    return(as.numeric(aggregated))
}

#' Expand hexbin results by metadata groups with feature expression
#'
#' Creates an expanded data frame where each bin appears once per metadata
#' group, with feature expression aggregated separately for each group.
#' Bins with zero cells for a particular group are excluded from that
#' group's facet.
#'
#' @param obj A SCUBA-supported single-cell object
#' @param cID Vector of cell-to-bin assignments from hexbin
#' @param hexbin_matrix Matrix with x, y coordinates and total cell counts
#' @param bin_ids Vector of unique bin IDs from hexbin object
#' @param feature_values Numeric vector of expression values (one per cell)
#' @param action Aggregation method: "mean", "median", or "sum"
#' @param metadata_var String naming the metadata variable
#' @param dr_name Name of dimension reduction (for cell ordering)
#' @param use_dims Dimensions used (for cell ordering)
#'
#' @return Data frame with columns: x, y, feature_value, and metadata_var
#' @import SCUBA
#' @importFrom stats setNames median
#' @keywords internal

.add_metadata_to_bins_with_feature <- function(obj, cID, hexbin_matrix, bin_ids, 
                                                feature_values, action, metadata_var, 
                                                dr_name, use_dims) {
    
    # Get cell IDs in the correct order
    dr_data <- fetch_reduction(obj, reduction = dr_name, dims = use_dims)
    cell_ids <- rownames(dr_data)
    
    # Fetch metadata for this variable
    metadata_df <- fetch_metadata(obj, vars = metadata_var)
    
    # Create cell_id -> metadata mapping
    cell_metadata <- setNames(metadata_df[[metadata_var]], rownames(metadata_df))
    
    # Match metadata to cells in the order they appear in the reduction
    cell_metadata_ordered <- cell_metadata[cell_ids]
    
    # Store original factor levels if this is a factor (for restoration later)
    is_factor <- is.factor(cell_metadata_ordered)
    if (is_factor) {
        original_levels <- levels(cell_metadata_ordered)
    }
    
    # Convert NAs to "NA" string for proper faceting (only if NAs actually exist)
    has_na <- any(is.na(cell_metadata_ordered))
    if (has_na) {
        if (is_factor) {
            levels(cell_metadata_ordered) <- c(levels(cell_metadata_ordered), "NA")
            original_levels <- levels(cell_metadata_ordered)  # Update stored levels
        }
        cell_metadata_ordered[is.na(cell_metadata_ordered)] <- "NA"
    }
    
    # Get unique bins that have at least one cell
    unique_bins <- unique(cID)
    
    # Create mapping from bin_id to row_index in hexbin_matrix
    bin_id_to_row <- setNames(seq_along(bin_ids), as.character(bin_ids))
    
    # Define aggregation function based on action
    agg_func <- switch(action,
        mean = function(x) mean(x, na.rm = TRUE),
        median = function(x) median(x, na.rm = TRUE),
        sum = function(x) sum(x, na.rm = TRUE)
    )
    
    # Create expanded data structure
    expanded_rows <- list()
    row_index <- 1
    
    for (bin_num in unique_bins) {
        # Find cells in this bin
        cells_in_bin <- which(cID == bin_num)
        
        # Get metadata values for cells in this bin
        meta_values <- cell_metadata_ordered[cells_in_bin]
        
        # Get feature values for cells in this bin
        feature_vals_in_bin <- feature_values[cells_in_bin]
        
        # Get unique groups in this bin
        unique_groups <- unique(meta_values)
        
        # Get row index for this bin ID
        row_idx <- bin_id_to_row[as.character(bin_num)]
        
        # Get bin coordinates using the row index
        bin_x <- hexbin_matrix[row_idx, "x"]
        bin_y <- hexbin_matrix[row_idx, "y"]
        
        # Create one row per group that has cells in this bin
        for (group in unique_groups) {
            # Get feature values for cells in this group within this bin
            group_mask <- meta_values == group
            group_feature_vals <- feature_vals_in_bin[group_mask]
            
            # Aggregate feature expression for this group
            agg_value <- agg_func(group_feature_vals)
            
            expanded_rows[[row_index]] <- data.frame(
                x = bin_x,
                y = bin_y,
                feature_value = agg_value,
                metadata_group = as.character(group),  # Convert to character to avoid factor level issues
                stringsAsFactors = FALSE
            )
            row_index <- row_index + 1
        }
    }
    
    # Combine all rows into a single data frame
    result <- do.call(rbind, expanded_rows)
    
    # Rename the metadata column to the actual variable name
    colnames(result)[colnames(result) == "metadata_group"] <- metadata_var
    
    # Restore factor type with original levels if input was a factor
    if (is_factor) {
        result[[metadata_var]] <- factor(result[[metadata_var]], levels = original_levels)
    }
    
    return(result)
}
