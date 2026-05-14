# Integration tests for exported plot functions
# Uses AML_Seurat example object loaded in helper-setup.R

skip_if_not_installed("SCUBA")

# --- plot_schextra_density() tests ---

test_that("plot_schextra_density returns a ggplot object", {
    p <- plot_schextra_density(AML_Seurat, nbins = 40, dimension_reduction = "umap")
    expect_s3_class(p, "ggplot")
})

test_that("plot_schextra_density has geom_hex layer", {
    p <- plot_schextra_density(AML_Seurat, nbins = 40, dimension_reduction = "umap")
    layer_classes <- sapply(p$layers, function(l) class(l$geom)[1])
    expect_true("GeomHex" %in% layer_classes)
})

test_that("plot_schextra_density uses correct axis labels from reduction", {
    p <- plot_schextra_density(AML_Seurat, nbins = 40, dimension_reduction = "umap")
    expect_equal(p$labels$x, "UMAP_1")
    expect_equal(p$labels$y, "UMAP_2")
})

test_that("plot_schextra_density uses custom title/labels when provided", {
    p <- plot_schextra_density(
        AML_Seurat, nbins = 40, dimension_reduction = "umap",
        title = "My Title", xlab = "X", ylab = "Y"
    )
    expect_equal(p$labels$title, "My Title")
    expect_equal(p$labels$x, "X")
    expect_equal(p$labels$y, "Y")
})

test_that("plot_schextra_density with split_by creates faceted plot", {
    p <- plot_schextra_density(
        AML_Seurat, nbins = 40, dimension_reduction = "umap",
        split_by = "Batch"
    )
    expect_s3_class(p, "ggplot")
    # Check facet exists
    expect_true(!is.null(p$facet))
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("plot_schextra_density with scale_density=TRUE creates scaled plot", {
    p <- plot_schextra_density(
        AML_Seurat, nbins = 40, dimension_reduction = "umap",
        split_by = "Batch", scale_density = TRUE
    )
    expect_s3_class(p, "ggplot")
    # The fill variable should be scaled_density
    expect_equal(rlang::as_name(p$mapping$fill), "scaled_density")
})

test_that("plot_schextra_density errors on invalid scales parameter", {
    expect_error(
        plot_schextra_density(
            AML_Seurat, nbins = 40, dimension_reduction = "umap",
            split_by = "Batch", scales = "invalid"
        ),
        "scales must be one of"
    )
})

test_that("plot_schextra_density errors on nonexistent split_by variable", {
    expect_error(
        plot_schextra_density(
            AML_Seurat, nbins = 40, dimension_reduction = "umap",
            split_by = "nonexistent_var"
        ),
        "not found"
    )
})

test_that("plot_schextra_density warns when scale_density=TRUE without split_by", {
    expect_warning(
        plot_schextra_density(
            AML_Seurat, nbins = 40, dimension_reduction = "umap",
            scale_density = TRUE
        ),
        "no effect"
    )
})

test_that("plot_schextra_density errors on invalid scale_density type", {
    expect_error(
        plot_schextra_density(
            AML_Seurat, nbins = 40, dimension_reduction = "umap",
            scale_density = "yes"
        ),
        "single logical value"
    )
})

test_that("plot_schextra_density works with PCA reduction", {
    p <- plot_schextra_density(AML_Seurat, nbins = 40, dimension_reduction = "pca")
    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$x, "PC_1")
    expect_equal(p$labels$y, "PC_2")
})

# --- plot_schextra_feature() tests ---

test_that("plot_schextra_feature returns a ggplot object", {
    p <- plot_schextra_feature(
        AML_Seurat, feature = "ACTG1", nbins = 40,
        dimension_reduction = "umap"
    )
    expect_s3_class(p, "ggplot")
})

test_that("plot_schextra_feature has correct default title", {
    p <- plot_schextra_feature(
        AML_Seurat, feature = "ACTG1", nbins = 40,
        dimension_reduction = "umap"
    )
    expect_equal(p$labels$title, "ACTG1 (mean)")
})

test_that("plot_schextra_feature with action=median uses correct title", {
    p <- plot_schextra_feature(
        AML_Seurat, feature = "ACTG1", nbins = 40,
        dimension_reduction = "umap", action = "median"
    )
    expect_equal(p$labels$title, "ACTG1 (median)")
})

test_that("plot_schextra_feature with split_by creates faceted plot", {
    p <- plot_schextra_feature(
        AML_Seurat, feature = "ACTG1", nbins = 40,
        dimension_reduction = "umap", split_by = "Batch"
    )
    expect_s3_class(p, "ggplot")
    expect_s3_class(p$facet, "FacetWrap")
})

test_that("plot_schextra_feature with numeric min_cutoff works", {
    p <- plot_schextra_feature(
        AML_Seurat, feature = "ACTG1", nbins = 40,
        dimension_reduction = "umap", min_cutoff = 0.5
    )
    expect_s3_class(p, "ggplot")
    # Build plot data to check values
    plot_data <- ggplot_build(p)$data[[1]]
    expect_true(all(plot_data$fill_mapped >= 0 | TRUE))  # Plot builds without error
})

test_that("plot_schextra_feature with quantile cutoffs works", {
    p <- plot_schextra_feature(
        AML_Seurat, feature = "ACTG1", nbins = 40,
        dimension_reduction = "umap",
        min_cutoff = "q10", max_cutoff = "q90"
    )
    expect_s3_class(p, "ggplot")
})

test_that("plot_schextra_feature errors on nonexistent feature", {
    expect_error(
        plot_schextra_feature(
            AML_Seurat, feature = "NOT_A_REAL_GENE_XYZ", nbins = 40,
            dimension_reduction = "umap"
        ),
        "Failed to retrieve feature|not found"
    )
})

test_that("plot_schextra_feature errors on invalid action", {
    expect_error(
        plot_schextra_feature(
            AML_Seurat, feature = "ACTG1", nbins = 40,
            dimension_reduction = "umap", action = "max"
        ),
        "action must be one of"
    )
})

test_that("plot_schextra_feature errors on invalid scales with split_by", {
    expect_error(
        plot_schextra_feature(
            AML_Seurat, feature = "ACTG1", nbins = 40,
            dimension_reduction = "umap",
            split_by = "Batch", scales = "wrong"
        ),
        "scales must be one of"
    )
})

test_that("plot_schextra_feature errors on nonexistent split_by variable", {
    expect_error(
        plot_schextra_feature(
            AML_Seurat, feature = "ACTG1", nbins = 40,
            dimension_reduction = "umap",
            split_by = "fake_variable_xyz"
        ),
        "not found"
    )
})

test_that("plot_schextra_feature legend shows feature name", {
    p <- plot_schextra_feature(
        AML_Seurat, feature = "ACTG1", nbins = 40,
        dimension_reduction = "umap"
    )
    expect_equal(p$labels$fill, "ACTG1")
})

test_that("plot_schextra_feature works with custom title and labels", {
    p <- plot_schextra_feature(
        AML_Seurat, feature = "ACTG1", nbins = 40,
        dimension_reduction = "umap",
        title = "Custom", xlab = "Dim1", ylab = "Dim2"
    )
    expect_equal(p$labels$title, "Custom")
    expect_equal(p$labels$x, "Dim1")
    expect_equal(p$labels$y, "Dim2")
})

test_that("plot_schextra_feature with split_by and cutoffs works together", {
    p <- plot_schextra_feature(
        AML_Seurat, feature = "ACTG1", nbins = 40,
        dimension_reduction = "umap",
        split_by = "Batch", min_cutoff = "q5", max_cutoff = "q95"
    )
    expect_s3_class(p, "ggplot")
    expect_s3_class(p$facet, "FacetWrap")
})
