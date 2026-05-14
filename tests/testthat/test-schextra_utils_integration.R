# Integration tests for internal utility functions requiring SCUBA
# Uses AML_Seurat example object loaded in helper-setup.R

skip_if_not_installed("SCUBA")

# --- .schextra_bin() tests ---

test_that(".schextra_bin returns correct list structure", {
    result <- test_bin_result
    expect_type(result, "list")
    expect_named(result, c("cID", "hexbin.matrix", "cell"))
})

test_that(".schextra_bin cID length matches cell count", {
    result <- test_bin_result
    n_cells <- ncol(AML_Seurat)
    expect_length(result$cID, n_cells)
})

test_that(".schextra_bin matrix has correct column names", {
    result <- test_bin_result
    expect_equal(colnames(result$hexbin.matrix), c("x", "y", "number_of_cells"))
})

test_that(".schextra_bin matrix rows match number of unique bins", {
    result <- test_bin_result
    expect_equal(nrow(result$hexbin.matrix), length(result$cell))
})

test_that(".schextra_bin all cID values map to valid bin IDs", {
    result <- test_bin_result
    expect_true(all(result$cID %in% result$cell))
})

test_that(".schextra_bin total cell count in matrix equals n_cells", {
    result <- test_bin_result
    expect_equal(sum(result$hexbin.matrix[, "number_of_cells"]), ncol(AML_Seurat))
})

test_that(".schextra_bin varying nbins changes bin count", {
    result_20 <- schextra:::.schextra_bin(AML_Seurat, nbins = 20, dr = "umap", use_dims = c(1, 2))
    result_80 <- schextra:::.schextra_bin(AML_Seurat, nbins = 80, dr = "umap", use_dims = c(1, 2))
    # More bins generally means more unique hexagons (not guaranteed but highly likely)
    expect_true(length(result_80$cell) >= length(result_20$cell))
})

test_that(".schextra_bin works with PCA reduction", {
    result <- schextra:::.schextra_bin(AML_Seurat, nbins = 30, dr = "pca", use_dims = c(1, 2))
    expect_type(result, "list")
    expect_equal(ncol(result$hexbin.matrix), 3)
    expect_length(result$cID, ncol(AML_Seurat))
})

# --- .add_metadata_to_bins() tests ---

test_that(".add_metadata_to_bins returns correct columns for character metadata", {
    result <- schextra:::.add_metadata_to_bins(
        AML_Seurat,
        test_bin_result$cID,
        test_bin_result$hexbin.matrix,
        test_bin_result$cell,
        "Batch",
        "umap",
        c(1, 2)
    )
    expect_s3_class(result, "data.frame")
    expect_true(all(c("x", "y", "number_of_cells", "Batch") %in% colnames(result)))
})

test_that(".add_metadata_to_bins per-bin group counts sum to original total", {
    result <- schextra:::.add_metadata_to_bins(
        AML_Seurat,
        test_bin_result$cID,
        test_bin_result$hexbin.matrix,
        test_bin_result$cell,
        "Batch",
        "umap",
        c(1, 2)
    )
    total_cells_expanded <- sum(result$number_of_cells)
    total_cells_original <- sum(test_bin_result$hexbin.matrix[, "number_of_cells"])
    expect_equal(total_cells_expanded, total_cells_original)
})

test_that(".add_metadata_to_bins preserves factor levels", {
    # condensed_cell_type is likely a character, so create a factor version
    AML_Seurat_factor <- AML_Seurat
    AML_Seurat_factor@meta.data$Batch_factor <- factor(
        AML_Seurat_factor@meta.data$Batch,
        levels = c("BM_200AB", "PBMC_200AB", "UNUSED_LEVEL")
    )
    
    bin_res <- schextra:::.schextra_bin(AML_Seurat_factor, nbins = 40, dr = "umap", use_dims = c(1, 2))
    
    result <- schextra:::.add_metadata_to_bins(
        AML_Seurat_factor,
        bin_res$cID,
        bin_res$hexbin.matrix,
        bin_res$cell,
        "Batch_factor",
        "umap",
        c(1, 2)
    )
    
    expect_true(is.factor(result$Batch_factor))
    # Original levels (minus unused) should be in the result levels
    expect_true("BM_200AB" %in% levels(result$Batch_factor))
    expect_true("PBMC_200AB" %in% levels(result$Batch_factor))
})

test_that(".add_metadata_to_bins handles NA metadata values", {
    # Inject NAs into a metadata column
    AML_Seurat_na <- AML_Seurat
    AML_Seurat_na@meta.data$Batch_na <- AML_Seurat_na@meta.data$Batch
    AML_Seurat_na@meta.data$Batch_na[1:10] <- NA
    
    bin_res <- schextra:::.schextra_bin(AML_Seurat_na, nbins = 40, dr = "umap", use_dims = c(1, 2))
    
    result <- schextra:::.add_metadata_to_bins(
        AML_Seurat_na,
        bin_res$cID,
        bin_res$hexbin.matrix,
        bin_res$cell,
        "Batch_na",
        "umap",
        c(1, 2)
    )
    
    # "NA" should appear as a group
    expect_true("NA" %in% result$Batch_na)
    # Total cells still matches
    expect_equal(sum(result$number_of_cells), ncol(AML_Seurat_na))
})

test_that(".add_metadata_to_bins only includes groups with cells in each bin", {
    result <- schextra:::.add_metadata_to_bins(
        AML_Seurat,
        test_bin_result$cID,
        test_bin_result$hexbin.matrix,
        test_bin_result$cell,
        "Batch",
        "umap",
        c(1, 2)
    )
    # All number_of_cells should be > 0 (no zero-count rows)
    expect_true(all(result$number_of_cells > 0))
})

# --- .add_metadata_to_bins_with_feature() tests ---

test_that(".add_metadata_to_bins_with_feature returns correct structure", {
    # Get feature values
    expr_data <- fetch_feature(AML_Seurat, features = "ACTG1")
    feature_values <- as.numeric(expr_data[["ACTG1"]])
    
    result <- schextra:::.add_metadata_to_bins_with_feature(
        AML_Seurat,
        test_bin_result$cID,
        test_bin_result$hexbin.matrix,
        test_bin_result$cell,
        feature_values,
        "mean",
        "Batch",
        "umap",
        c(1, 2)
    )
    
    expect_s3_class(result, "data.frame")
    expect_true(all(c("x", "y", "feature_value", "Batch") %in% colnames(result)))
})

test_that(".add_metadata_to_bins_with_feature aggregates with sum correctly", {
    expr_data <- fetch_feature(AML_Seurat, features = "ACTG1")
    feature_values <- as.numeric(expr_data[["ACTG1"]])
    
    result <- schextra:::.add_metadata_to_bins_with_feature(
        AML_Seurat,
        test_bin_result$cID,
        test_bin_result$hexbin.matrix,
        test_bin_result$cell,
        feature_values,
        "sum",
        "Batch",
        "umap",
        c(1, 2)
    )
    
    # Sum of all per-group feature sums should equal total sum
    # (since each cell appears in exactly one group)
    expect_equal(sum(result$feature_value), sum(feature_values), tolerance = 1e-10)
})

test_that(".add_metadata_to_bins_with_feature handles median action", {
    expr_data <- fetch_feature(AML_Seurat, features = "ACTG1")
    feature_values <- as.numeric(expr_data[["ACTG1"]])
    
    result <- schextra:::.add_metadata_to_bins_with_feature(
        AML_Seurat,
        test_bin_result$cID,
        test_bin_result$hexbin.matrix,
        test_bin_result$cell,
        feature_values,
        "median",
        "Batch",
        "umap",
        c(1, 2)
    )
    
    # All feature_values should be finite (no NaN from empty groups)
    expect_true(all(is.finite(result$feature_value)))
})

test_that(".add_metadata_to_bins_with_feature excludes zero-cell groups", {
    expr_data <- fetch_feature(AML_Seurat, features = "ACTG1")
    feature_values <- as.numeric(expr_data[["ACTG1"]])
    
    result <- schextra:::.add_metadata_to_bins_with_feature(
        AML_Seurat,
        test_bin_result$cID,
        test_bin_result$hexbin.matrix,
        test_bin_result$cell,
        feature_values,
        "mean",
        "Batch",
        "umap",
        c(1, 2)
    )
    
    # Every row should represent a real group with cells
    # (no NaN from mean of empty vector, though na.rm handles it)
    expect_true(nrow(result) > 0)
    expect_true(all(!is.nan(result$feature_value)))
})
