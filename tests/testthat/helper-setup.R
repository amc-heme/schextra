# Helper file auto-loaded by testthat before tests run
# Provides shared fixtures for tests that require SCUBA

# Load the AML_Seurat example object if SCUBA is available
if (requireNamespace("SCUBA", quietly = TRUE)) {
    library(SCUBA)
    # Load into a temporary environment then assign to this scope
    tmp_env <- new.env(parent = emptyenv())
    data(AML_Seurat, package = "SCUBA", envir = tmp_env)
    AML_Seurat <- tmp_env$AML_Seurat
    
    # Pre-compute a binning result for reuse across tests
    # Using smaller nbins for speed in tests
    test_bin_result <- schextra:::.schextra_bin(
        AML_Seurat, nbins = 40, dr = "umap", use_dims = c(1, 2)
    )
}
