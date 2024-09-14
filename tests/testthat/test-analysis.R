# Test standard.p.value function
test_that("standard.p.value returns correct text for small values", {
  expect_equal(standard.p.value(0.0001, 3), "<0.001")
  expect_equal(standard.p.value(0.00001, 4), "<1e-04")
})

test_that("standard.p.value rounds values correctly", {
  expect_equal(standard.p.value(0.123456, 3), 0.123)
  expect_equal(standard.p.value(0.123456, 4), 0.1235)
})

test_that("standard.p.value handles boundary cases", {
  expect_equal(standard.p.value(0.001, 3), 0.001)
  expect_equal(standard.p.value(0.00099, 3), "<0.001")
})
