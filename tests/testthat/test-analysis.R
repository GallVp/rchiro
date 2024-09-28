# Test describe_p_value function
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
