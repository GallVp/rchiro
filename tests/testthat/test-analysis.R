# is.nan.data.frame ---------------------------------------------------------------------------------------------------
test_that("is.nan.data.frame works correctly", {
  # Create a data frame with numeric values, including NaN
  df <- data.frame(
    a = c(1, 2, NaN, 4),
    b = c(NaN, 5, 6, NaN)
  )

  # Expected output: a logical matrix indicating where NaNs are in the data frame
  expected_output <- data.frame(
    a = c(FALSE, FALSE, TRUE, FALSE),
    b = c(TRUE, FALSE, FALSE, TRUE)
  )

  # Test if the function works as expected
  result <- is.nan.data.frame(df)

  expect_equal(result, expected_output)
})

test_that("is.nan.data.frame handles data frames without NaN", {
  # Create a data frame without NaN values
  df <- data.frame(
    a = c(1, 2, 3, 4),
    b = c(5, 6, 7, 8)
  )

  # Expected output: all FALSE since there are no NaN values
  expected_output <- data.frame(
    a = c(FALSE, FALSE, FALSE, FALSE),
    b = c(FALSE, FALSE, FALSE, FALSE)
  )

  # Test if the function works as expected
  result <- is.nan.data.frame(df)

  expect_equal(result, expected_output)
})

test_that("is.nan.data.frame handles empty data frames", {
  # Create an empty data frame
  df <- data.frame()

  # Expected output: an empty data frame
  expected_output <- data.frame()

  # Test if the function works as expected
  result <- is.nan.data.frame(df)

  expect_equal(result, expected_output)
})

# describe_p_value ----------------------------------------------------------------------------------------------------
test_that("describe_p_value returns correct text for small values", {
  expect_equal(describe_p_value(0.0001, 3), "<0.001")
  expect_equal(describe_p_value(0.00001, 4), "<1e-04")
})

test_that("describe_p_value rounds values correctly", {
  expect_equal(describe_p_value(0.123456, 3), 0.123)
  expect_equal(describe_p_value(0.123456, 4), 0.1235)
})

test_that("describe_p_value handles boundary cases", {
  expect_equal(describe_p_value(0.001, 3), 0.001)
  expect_equal(describe_p_value(0.00099, 3), "<0.001")
})

# names_of_df_except --------------------------------------------------------------------------------------------------
test_that("names_of_df_except works correctly with a single column exclusion", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)

  # Exclude column 'b'
  result <- names_of_df_except(df, "b")

  # Expected result
  expected_result <- c("a", "c")

  expect_equal(result, expected_result)
})

test_that("names_of_df_except works with multiple column exclusions", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9, d = 10:12)

  # Exclude columns 'b' and 'c'
  result <- names_of_df_except(df, c("b", "c"))

  # Expected result
  expected_result <- c("a", "d")

  expect_equal(result, expected_result)
})

test_that("names_of_df_except returns all names when no exclusions", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)

  # No exclusions
  result <- names_of_df_except(df, character(0))

  # Expected result: all column names
  expected_result <- c("a", "b", "c")

  expect_equal(result, expected_result)
})

test_that("names_of_df_except returns empty vector when all columns are excluded", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)

  # Exclude all columns
  result <- names_of_df_except(df, c("a", "b", "c"))

  # Expected result: empty character vector
  expected_result <- character(0)

  expect_equal(result, expected_result)
})

# columns_of_df_except ------------------------------------------------------------------------------------------------
test_that("columns_of_df_except works with a single column exclusion", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)

  # Exclude column 'b'
  result <- columns_of_df_except(df, "b")

  # Expected result: data frame without column 'b'
  expected_result <- data.frame(a = 1:3, c = 7:9)

  expect_equal(result, expected_result)
})

test_that("columns_of_df_except works with multiple column exclusions", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9, d = 10:12)

  # Exclude columns 'b' and 'c'
  result <- columns_of_df_except(df, c("b", "c"))

  # Expected result: data frame without columns 'b' and 'c'
  expected_result <- data.frame(a = 1:3, d = 10:12)

  expect_equal(result, expected_result)
})

test_that("columns_of_df_except returns full data frame when no exclusions", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)

  # No exclusions
  result <- columns_of_df_except(df, character(0))

  # Expected result: return full data frame
  expected_result <- df

  expect_equal(result, expected_result)
})

test_that("columns_of_df_except returns empty data frame when all columns are excluded", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)

  # Exclude all columns
  result <- columns_of_df_except(df, c("a", "b", "c"))

  # Expected result: an empty data frame with correct number of rows
  expected_result <- df[, character(0)]

  expect_equal(result, expected_result)
})

# separate_baseline_df ------------------------------------------------------------------------------------------------
test_that("separate_baseline_df correctly separates baseline and merges with post-baseline", {
  # Simulate a simple dataset
  df <- data.frame(
    PartId = rep(1:3, each = 2),
    Group = rep(c("A", "B", "A"), each = 2),
    Time = c("Baseline", "Follow-up", "Baseline", "Follow-up", "Baseline", "Follow-up"),
    Outcome = c(10, 15, 20, 25, 30, 35)
  )

  # Call the function
  result <- separate_baseline_df(df)

  # Expected result
  expected_result <- data.frame(
    PartId = rep(1:3, each = 1),
    Group = c("A", "B", "A"),
    Time = c("Follow-up", "Follow-up", "Follow-up"),
    Outcome = c(15, 25, 35),
    OutcomePre = c(10, 20, 30)
  )

  expect_equal(result, expected_result)
})

test_that("separate_baseline_df handles different time variable names and baseline levels", {
  # Simulate a dataset with different variable names
  df <- data.frame(
    PartId = rep(1:3, each = 2),
    Group = rep(c("A", "B", "A"), each = 2),
    TimePoint = c("Pre", "Post", "Pre", "Post", "Pre", "Post"),
    Score = c(10, 15, 20, 25, 30, 35)
  )

  # Call the function with custom parameters
  result <- separate_baseline_df(df, time_var = "TimePoint", baseline_lvl = "Pre", var_name = "Score")

  # Expected result
  expected_result <- data.frame(
    PartId = rep(1:3, each = 1),
    Group = c("A", "B", "A"),
    TimePoint = c("Post", "Post", "Post"),
    Score = c(15, 25, 35),
    ScorePre = c(10, 20, 30)
  )

  expect_equal(result, expected_result)
})

test_that("separate_baseline_df handles cases with multiple post-baseline observations", {
  # Simulate a dataset with multiple post-baseline observations
  df <- data.frame(
    PartId = rep(1:2, each = 3),
    Group = rep(c("A", "B"), each = 3),
    Time = c("Baseline", "Follow-up1", "Follow-up2", "Baseline", "Follow-up1", "Follow-up2"),
    Outcome = c(10, 15, 20, 30, 35, 40)
  )

  # Call the function
  result <- separate_baseline_df(df)

  # Expected result
  expected_result <- data.frame(
    PartId = rep(1:2, each = 2),
    Group = c("A", "A", "B", "B"),
    Time = c("Follow-up1", "Follow-up2", "Follow-up1", "Follow-up2"),
    Outcome = c(15, 20, 35, 40),
    OutcomePre = c(10, 10, 30, 30)
  )

  expect_equal(result, expected_result)
})
