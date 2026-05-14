# Multi-backend integration tests
# Verifies schextra works across all SCUBA-supported object types
# Uses DRY helper to run core test battery per backend

skip_if_not_installed("SCUBA")

# DRY helper: runs core test battery for any backend
# Parameters:
#   obj - the single-cell object
#   label - string prefix for test names (e.g., "SCE")
#   dr - reduction name (e.g., "umap", "UMAP", "X_umap")
#   test_feature - logical, whether to test feature functions
run_backend_tests <- function(obj, label, dr, test_feature = TRUE) {

    # -- Binning tests --
    test_that(paste(label, "- .schextra_bin returns valid structure"), {
        res <- schextra:::.schextra_bin(obj, nbins = 40, dr = dr, use_dims = c(1, 2))
        expect_type(res, "list")
        expect_named(res, c("cID", "hexbin.matrix", "cell"))
        expect_equal(colnames(res$hexbin.matrix), c("x", "y", "number_of_cells"))
    })

    test_that(paste(label, "- .schextra_bin cID length matches cell count"), {
        res <- schextra:::.schextra_bin(obj, nbins = 40, dr = dr, use_dims = c(1, 2))
        expect_length(res$cID, 250)
    })

    # -- Density plot tests --
    test_that(paste(label, "- plot_schextra_density returns ggplot"), {
        p <- plot_schextra_density(obj, nbins = 40, dimension_reduction = dr)
        expect_s3_class(p, "ggplot")
    })

    test_that(paste(label, "- plot_schextra_density with split_by works"), {
        p <- plot_schextra_density(obj, nbins = 40, dimension_reduction = dr, split_by = "Batch")
        expect_s3_class(p, "ggplot")
        expect_s3_class(p$facet, "FacetWrap")
    })

    # -- Metadata bin expansion test --
    test_that(paste(label, "- .add_metadata_to_bins returns correct structure"), {
        bin_res <- schextra:::.schextra_bin(obj, nbins = 40, dr = dr, use_dims = c(1, 2))
        result <- schextra:::.add_metadata_to_bins(
            obj, bin_res$cID, bin_res$hexbin.matrix, bin_res$cell,
            "Batch", dr, c(1, 2)
        )
        expect_s3_class(result, "data.frame")
        expect_true(all(c("x", "y", "number_of_cells", "Batch") %in% colnames(result)))
        expect_true(all(result$number_of_cells > 0))
        expect_equal(sum(result$number_of_cells), 250)
    })

    # -- Feature tests (only when fetch_feature works for this backend) --
    if (test_feature) {
        test_that(paste(label, "- plot_schextra_feature returns ggplot"), {
            p <- plot_schextra_feature(obj, feature = "ACTG1", nbins = 40,
                                       dimension_reduction = dr)
            expect_s3_class(p, "ggplot")
        })

        test_that(paste(label, "- plot_schextra_feature with split_by works"), {
            p <- plot_schextra_feature(obj, feature = "ACTG1", nbins = 40,
                                       dimension_reduction = dr, split_by = "Batch")
            expect_s3_class(p, "ggplot")
            expect_s3_class(p$facet, "FacetWrap")
        })

        test_that(paste(label, "- .add_metadata_to_bins_with_feature returns correct structure"), {
            bin_res <- schextra:::.schextra_bin(obj, nbins = 40, dr = dr, use_dims = c(1, 2))
            expr_data <- SCUBA::fetch_feature(obj, features = "ACTG1")
            feature_values <- as.numeric(expr_data[["ACTG1"]])

            result <- schextra:::.add_metadata_to_bins_with_feature(
                obj, bin_res$cID, bin_res$hexbin.matrix, bin_res$cell,
                feature_values, "mean", "Batch", dr, c(1, 2)
            )
            expect_s3_class(result, "data.frame")
            expect_true(all(c("x", "y", "feature_value", "Batch") %in% colnames(result)))
            expect_true(all(!is.nan(result$feature_value)))
        })
    }
}

# --- SingleCellExperiment backend ---
# Reduction name: "UMAP" (uppercase)
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    sce_obj <- SCUBA::AML_SCE()
    run_backend_tests(sce_obj, "SCE", dr = "UMAP", test_feature = TRUE)
}

# --- AnnData backend (in-memory h5ad) ---
# Reduction name: "X_umap" (anndata obsm key convention)
if (requireNamespace("anndata", quietly = TRUE) &&
    requireNamespace("reticulate", quietly = TRUE) &&
    reticulate::py_available()) {
    h5ad_obj <- SCUBA::AML_h5ad()
    run_backend_tests(h5ad_obj, "AnnData", dr = "X_umap", test_feature = TRUE)
}

# --- BPCells-backed Seurat ---
# Reduction name: "umap" (lowercase, same as regular Seurat)
# NOTE: test_feature = FALSE because the bundled AML_BPCells object has a
# hardcoded backing store path that doesn't exist outside the package
# developer's machine, causing fetch_feature() to fail.
if (requireNamespace("BPCells", quietly = TRUE)) {
    bp_env <- new.env(parent = emptyenv())
    data(AML_BPCells, package = "SCUBA", envir = bp_env)
    bp_obj <- bp_env$AML_BPCells
    run_backend_tests(bp_obj, "BPCells", dr = "umap", test_feature = FALSE)
}
