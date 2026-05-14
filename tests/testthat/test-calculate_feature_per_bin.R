# Tests for .calculate_feature_per_bin()
# Pure function -- no SCUBA dependency

test_that(".calculate_feature_per_bin computes mean correctly", {
    # 3 cells in bin 1 (values 2,4,6 -> mean=4), 2 cells in bin 2 (values 10,20 -> mean=15)
    feature_values <- c(2, 4, 6, 10, 20)
    cID <- c(1, 1, 1, 2, 2)
    result <- schextra:::.calculate_feature_per_bin(feature_values, cID, "mean")
    expect_equal(result, c(4, 15))
})

test_that(".calculate_feature_per_bin computes median correctly", {
    feature_values <- c(1, 5, 9, 10, 30)
    cID <- c(1, 1, 1, 2, 2)
    result <- schextra:::.calculate_feature_per_bin(feature_values, cID, "median")
    expect_equal(result, c(5, 20))  # median(1,5,9)=5, median(10,30)=20
})

test_that(".calculate_feature_per_bin computes sum correctly", {
    feature_values <- c(2, 3, 5, 10, 20)
    cID <- c(1, 1, 1, 2, 2)
    result <- schextra:::.calculate_feature_per_bin(feature_values, cID, "sum")
    expect_equal(result, c(10, 30))  # sum(2,3,5)=10, sum(10,20)=30
})

test_that(".calculate_feature_per_bin errors on invalid action", {
    expect_error(
        schextra:::.calculate_feature_per_bin(c(1, 2), c(1, 1), "max"),
        "action must be one of"
    )
})

test_that(".calculate_feature_per_bin errors on non-numeric feature_values", {
    expect_error(
        schextra:::.calculate_feature_per_bin(c("a", "b"), c(1, 1), "mean"),
        "feature_values must be numeric"
    )
})

test_that(".calculate_feature_per_bin handles NAs with na.rm", {
    feature_values <- c(2, NA, 6, NA, 20)
    cID <- c(1, 1, 1, 2, 2)
    result <- schextra:::.calculate_feature_per_bin(feature_values, cID, "mean")
    expect_equal(result, c(4, 20))  # mean(2,NA,6)=4, mean(NA,20)=20
})

test_that(".calculate_feature_per_bin handles single cell per bin", {
    feature_values <- c(5, 10, 15)
    cID <- c(1, 2, 3)
    result <- schextra:::.calculate_feature_per_bin(feature_values, cID, "mean")
    expect_equal(result, c(5, 10, 15))
})

test_that(".calculate_feature_per_bin output length equals unique bins", {
    feature_values <- c(1, 2, 3, 4, 5, 6, 7, 8)
    cID <- c(10, 10, 20, 20, 20, 30, 30, 30)
    result <- schextra:::.calculate_feature_per_bin(feature_values, cID, "sum")
    expect_length(result, 3)
})

test_that(".calculate_feature_per_bin handles all-zero values", {
    feature_values <- c(0, 0, 0, 0)
    cID <- c(1, 1, 2, 2)
    result <- schextra:::.calculate_feature_per_bin(feature_values, cID, "mean")
    expect_equal(result, c(0, 0))
})
