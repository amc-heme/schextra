# Tests for .apply_cutoffs()
# Pure function -- no SCUBA dependency

test_that(".apply_cutoffs returns unchanged values with NULL cutoffs", {
    vals <- c(1, 5, 10, 15, 20)
    result <- schextra:::.apply_cutoffs(vals)
    expect_equal(result, vals)
})

test_that(".apply_cutoffs clamps low values with numeric min_cutoff", {
    vals <- c(1, 5, 10, 15, 20)
    result <- schextra:::.apply_cutoffs(vals, min_cutoff = 5)
    expect_equal(result, c(5, 5, 10, 15, 20))
})

test_that(".apply_cutoffs clamps high values with numeric max_cutoff", {
    vals <- c(1, 5, 10, 15, 20)
    result <- schextra:::.apply_cutoffs(vals, max_cutoff = 15)
    expect_equal(result, c(1, 5, 10, 15, 15))
})

test_that(".apply_cutoffs applies both numeric cutoffs", {
    vals <- c(1, 5, 10, 15, 20)
    result <- schextra:::.apply_cutoffs(vals, min_cutoff = 5, max_cutoff = 15)
    expect_equal(result, c(5, 5, 10, 15, 15))
})

test_that(".apply_cutoffs handles quantile string for min_cutoff", {
    vals <- 1:100
    result <- schextra:::.apply_cutoffs(vals, min_cutoff = "q10")
    # 10th percentile of 1:100 = 10.9
    cutoff <- quantile(vals, probs = 0.10)
    expect_true(all(result >= cutoff))
    expect_equal(min(result), as.numeric(cutoff))
    # Values above cutoff unchanged
    expect_equal(result[result > cutoff], vals[vals > cutoff])
})

test_that(".apply_cutoffs handles quantile string for max_cutoff", {
    vals <- 1:100
    result <- schextra:::.apply_cutoffs(vals, max_cutoff = "q90")
    cutoff <- quantile(vals, probs = 0.90)
    expect_true(all(result <= cutoff))
    expect_equal(max(result), as.numeric(cutoff))
    # Values below cutoff unchanged
    expect_equal(result[result < cutoff], vals[vals < cutoff])
})

test_that(".apply_cutoffs handles both quantile cutoffs", {
    vals <- 1:100
    result <- schextra:::.apply_cutoffs(vals, min_cutoff = "q5", max_cutoff = "q95")
    min_c <- quantile(vals, probs = 0.05)
    max_c <- quantile(vals, probs = 0.95)
    expect_true(all(result >= min_c))
    expect_true(all(result <= max_c))
})

test_that(".apply_cutoffs errors on out-of-range quantile", {
    vals <- 1:10
    expect_error(
        schextra:::.apply_cutoffs(vals, min_cutoff = "q101"),
        "between 0 and 100"
    )
    expect_error(
        schextra:::.apply_cutoffs(vals, max_cutoff = "q101"),
        "between 0 and 100"
    )
})

test_that(".apply_cutoffs errors on invalid string format", {
    vals <- 1:10
    expect_error(
        schextra:::.apply_cutoffs(vals, min_cutoff = "abc"),
        "must be NULL, a numeric value, or a quantile string"
    )
    expect_error(
        schextra:::.apply_cutoffs(vals, max_cutoff = "foo"),
        "must be NULL, a numeric value, or a quantile string"
    )
})

test_that(".apply_cutoffs errors on invalid type (logical)", {
    vals <- 1:10
    expect_error(
        schextra:::.apply_cutoffs(vals, min_cutoff = TRUE),
        "must be NULL, a numeric value, or a quantile string"
    )
    expect_error(
        schextra:::.apply_cutoffs(vals, max_cutoff = FALSE),
        "must be NULL, a numeric value, or a quantile string"
    )
})

test_that(".apply_cutoffs handles NAs in input values", {
    vals <- c(1, NA, 5, NA, 10)
    result <- schextra:::.apply_cutoffs(vals, min_cutoff = 3, max_cutoff = 8)
    # pmax/pmin with na.rm=TRUE treats NAs as -Inf/Inf respectively,
    # so NAs get replaced by the cutoff value
    expect_equal(result, c(3, 3, 5, 3, 8))
})

test_that(".apply_cutoffs handles min > max (clamps to narrow range)", {
    vals <- c(1, 5, 10, 15, 20)
    # min_cutoff=12, max_cutoff=8: after pmax(vals,12) -> c(12,12,12,15,20)
    # then pmin(result,8) -> c(8,8,8,8,8)
    result <- schextra:::.apply_cutoffs(vals, min_cutoff = 12, max_cutoff = 8)
    expect_equal(result, c(8, 8, 8, 8, 8))
})

test_that(".apply_cutoffs q0 and q100 are valid edge cases", {
    vals <- c(2, 5, 8, 12, 20)
    # q0 = minimum value, q100 = maximum value
    result_min <- schextra:::.apply_cutoffs(vals, min_cutoff = "q0")
    expect_equal(result_min, vals)  # q0 = min(vals) = 2, no change
    
    result_max <- schextra:::.apply_cutoffs(vals, max_cutoff = "q100")
    expect_equal(result_max, vals)  # q100 = max(vals) = 20, no change
})
