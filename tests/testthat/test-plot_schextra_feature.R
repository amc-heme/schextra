# Tests for plot_schextra_feature(), focused on the metadata fallback fix
# (numeric metadata variables such as nCount_RNA) while preserving the
# BPCells direct-matrix fast path.
#
# Uses the small example objects shipped with SCUBA:
#   - AML_Seurat:  in-memory Seurat object (RNA default assay, AB second assay)
#   - AML_BPCells: BPCells-backed Seurat object (RNA/AB IterableMatrix layers)

skip_if_no_scuba <- function() {
    testthat::skip_if_not_installed("SCUBA")
    testthat::skip_if_not_installed("SeuratObject")
}

# Load AML_Seurat from SCUBA, skipping if unavailable.
load_aml_seurat <- function() {
    skip_if_no_scuba()
    e <- new.env()
    loaded <- tryCatch({
        utils::data("AML_Seurat", package = "SCUBA", envir = e)
        TRUE
    }, error = function(err) FALSE)
    if (!loaded || !exists("AML_Seurat", envir = e)) {
        testthat::skip("AML_Seurat example object not available")
    }
    get("AML_Seurat", envir = e)
}

# Load AML_BPCells from SCUBA and repair the on-disk matrix directory paths.
# The shipped object stores absolute paths from the package author's machine;
# set_bpcells_dir() repoints each assay/layer to the installed extdata copy so
# the BPCells direct-matrix path can actually be exercised.
load_aml_bpcells <- function() {
    skip_if_no_scuba()
    testthat::skip_if_not_installed("BPCells")
    e <- new.env()
    loaded <- tryCatch({
        utils::data("AML_BPCells", package = "SCUBA", envir = e)
        TRUE
    }, error = function(err) FALSE)
    if (!loaded || !exists("AML_BPCells", envir = e)) {
        testthat::skip("AML_BPCells example object not available")
    }
    obj <- get("AML_BPCells", envir = e)

    base <- system.file("extdata", "AML_BPCells", package = "SCUBA")
    if (!nzchar(base) || !dir.exists(base)) {
        testthat::skip("AML_BPCells extdata matrices not available")
    }
    for (a in c("RNA", "AB")) {
        for (l in c("counts", "data")) {
            d <- file.path(base, paste0(a, "_", l))
            if (dir.exists(d)) {
                obj <- SCUBA::set_bpcells_dir(
                    obj, assay = a, layer = l, dirname = d
                )
            }
        }
    }
    obj
}


test_that("case 1: plots a gene from the default assay (unkeyed name)", {
    obj <- load_aml_seurat()

    p <- plot_schextra_feature(
        obj,
        feature = "ACTG1",
        dimension_reduction = "umap"
    )

    expect_s3_class(p, "ggplot")
    vals <- p$data$feature_value
    expect_true(length(vals) > 0)
    expect_true(all(is.finite(vals)))
})


test_that("case 2: plots a keyed feature from a non-default assay", {
    obj <- load_aml_seurat()

    # "ab_CD10-AB" is the keyed name (AB assay key "ab_" + feature "CD10-AB").
    p <- plot_schextra_feature(
        obj,
        feature = "ab_CD10-AB",
        assay = "AB",
        dimension_reduction = "umap"
    )

    expect_s3_class(p, "ggplot")
    expect_true(length(p$data$feature_value) > 0)
    expect_true(all(is.finite(p$data$feature_value)))
})


test_that("case 3: plots a numeric metadata variable (nCount_RNA)", {
    obj <- load_aml_seurat()

    # Previously failing case: fetch_feature cannot resolve metadata names.
    expect_no_error(
        p <- plot_schextra_feature(
            obj,
            feature = "nCount_RNA",
            dimension_reduction = "umap"
        )
    )

    expect_s3_class(p, "ggplot")
    expect_true(length(p$data$feature_value) > 0)
    expect_true(all(is.finite(p$data$feature_value)))
})


test_that("case 4: numeric metadata variable works with split_by", {
    obj <- load_aml_seurat()

    expect_no_error(
        p <- plot_schextra_feature(
            obj,
            feature = "nCount_RNA",
            split_by = "condensed_cell_type",
            dimension_reduction = "umap"
        )
    )

    expect_s3_class(p, "ggplot")
    # split_by exercises .add_metadata_to_bins_with_feature(); the resulting
    # data frame must carry the split variable used for faceting.
    expect_true("condensed_cell_type" %in% colnames(p$data))
    expect_true(length(p$data$feature_value) > 0)
    expect_true(all(is.finite(p$data$feature_value)))
})


test_that("case 5: BPCells gene uses the direct-matrix path; metadata works", {
    obj <- load_aml_bpcells()

    # The fast path is taken when the assay/layer is BPCells-backed. If it is
    # not, fetch_feature emits "Falling back to fetch_data." -- assert both that
    # the assay is BPCells-backed and that no fallback message is emitted.
    expect_true(SCUBA:::is_bpcells(obj, assay = "RNA", layer = "data"))

    msgs <- character(0)
    withCallingHandlers(
        p <- plot_schextra_feature(
            obj,
            feature = "ACTG1",
            assay = "RNA",
            layer = "data",
            dimension_reduction = "umap"
        ),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )

    expect_s3_class(p, "ggplot")
    expect_false(any(grepl("Falling back to fetch_data", msgs, fixed = TRUE)))
    expect_true(length(p$data$feature_value) > 0)
    expect_true(all(is.finite(p$data$feature_value)))

    # Compare against direct fetch_feature retrieval to confirm correctness of
    # the values flowing through the fast path.
    expr <- SCUBA::fetch_feature(
        obj, features = "ACTG1", assay = "RNA", layer = "data"
    )
    expect_true(all(is.finite(as.numeric(expr[["ACTG1"]]))))

    # Metadata fallback must also work on the BPCells-backed object.
    expect_no_error(
        p_meta <- plot_schextra_feature(
            obj,
            feature = "nCount_RNA",
            dimension_reduction = "umap"
        )
    )
    expect_s3_class(p_meta, "ggplot")
    expect_true(all(is.finite(p_meta$data$feature_value)))
})


test_that("case 6: a nonexistent name errors with the existing message", {
    obj <- load_aml_seurat()

    expect_error(
        plot_schextra_feature(
            obj,
            feature = "NOT_A_REAL_FEATURE_XYZ",
            dimension_reduction = "umap"
        ),
        regexp = "Failed to retrieve feature|not found in the specified assay|not found in assay"
    )
})


# ---------------------------------------------------------------------------
# Resolution semantics under a gene/metadata name collision.
#
# We inject a metadata column named exactly like a real gene ("ACTG1") and
# give it distinct, easily distinguishable values (1..N, one per cell). Real
# ACTG1 expression in this object is small (~0-6), so a per-bin mean cleanly
# separates the two sources:
#   - metadata source -> bin means are large (tens to hundreds)
#   - gene source     -> bin means are small (< 10)
# ---------------------------------------------------------------------------

make_collision_obj <- function() {
    obj <- load_aml_seurat()
    testthat::skip_if_not("ACTG1" %in% rownames(obj[["RNA"]]))
    # Inject a metadata column colliding with the gene name.
    obj@meta.data[["ACTG1"]] <- seq_len(ncol(obj))
    obj
}

# Largest per-bin feature value in the plot; used to discriminate source.
max_bin_value <- function(p) {
    max(p$data$feature_value, na.rm = TRUE)
}

test_that("collision: assay supplied returns the GENE (assay-only)", {
    obj <- make_collision_obj()

    p <- plot_schextra_feature(
        obj,
        feature = "ACTG1",
        assay = "RNA",
        dimension_reduction = "umap"
    )

    expect_s3_class(p, "ggplot")
    # Gene expression is small; injected metadata (1..N) would be large.
    expect_lt(max_bin_value(p), 50)
})

test_that("collision: no assay returns METADATA (metadata-first default)", {
    obj <- make_collision_obj()

    # SCUBA emits a warning when a bare name collides with metadata and the
    # default assay, telling the user to use the key for the gene. That
    # warning is expected here (we are deliberately resolving the metadata).
    suppressWarnings(
        p <- plot_schextra_feature(
            obj,
            feature = "ACTG1",
            dimension_reduction = "umap"
        )
    )

    expect_s3_class(p, "ggplot")
    # Metadata values are 1..N (N = number of cells), so per-bin means are
    # far larger than any real ACTG1 expression.
    expect_gt(max_bin_value(p), 50)
})

test_that("collision: keyed name returns the GENE with or without assay", {
    obj <- make_collision_obj()

    p_no_assay <- plot_schextra_feature(
        obj,
        feature = "rna_ACTG1",
        dimension_reduction = "umap"
    )
    p_assay <- plot_schextra_feature(
        obj,
        feature = "rna_ACTG1",
        assay = "RNA",
        dimension_reduction = "umap"
    )

    expect_s3_class(p_no_assay, "ggplot")
    expect_s3_class(p_assay, "ggplot")
    expect_lt(max_bin_value(p_no_assay), 50)
    expect_lt(max_bin_value(p_assay), 50)
})

test_that("collision: assay supplied + metadata-only name errors (strict)", {
    obj <- make_collision_obj()

    # nCount_RNA is metadata only (not an RNA gene); with assay = "RNA" the
    # search is restricted to the assay and must error rather than resolve
    # from metadata.
    expect_error(
        plot_schextra_feature(
            obj,
            feature = "nCount_RNA",
            assay = "RNA",
            dimension_reduction = "umap"
        ),
        regexp = "not found in assay"
    )
})
