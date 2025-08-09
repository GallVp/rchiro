# Tests for describe_contrasts and describe_emmeans

test_that("describe_contrasts returns expected columns and formatting", {
  library(emmeans)
  library(dplyr)
  model <- lm(Sepal.Length ~ Sepal.Width + Species, data = iris)
  contrast_table <- emmeans(model, pairwise ~ Species)$contrasts
  result <- rchiro::describe_contrasts(contrast_table)
  expect_true(is.data.frame(result))
  expect_true(any(grepl("Difference±SE", names(result))))
  expect_true(any(grepl("p-value", names(result))))
  expect_true(any(grepl("Contrast", names(result))))
  expect_true(all(grepl("±", result[[grep("Difference±SE", names(result))]])))
  expect_true(all(grepl("[0-9]", result[[grep("p-value", names(result))]])))
})

test_that("describe_contrasts handles asymptotic and regular CIs", {
  library(emmeans)
  model <- lm(Sepal.Length ~ Sepal.Width + Species, data = iris)
  contrast_table <- emmeans(model, pairwise ~ Species)$contrasts
  # Should work for regular CIs
  expect_silent(rchiro::describe_contrasts(contrast_table))
  # Should work for summary with infer=TRUE (asymptotic)
  expect_silent(rchiro::describe_contrasts(summary(contrast_table, infer = c(TRUE, TRUE))))
})

test_that("describe_emmeans returns expected columns and formatting", {
  library(emmeans)
  model <- lm(Sepal.Length ~ Sepal.Width + Species, data = iris)
  emm_table <- emmeans(model, ~Species)
  result <- rchiro::describe_emmeans(emm_table)
  expect_true(is.data.frame(result))
  expect_true(any(grepl("Estimate±SE", names(result))))
  expect_true(any(grepl("p-value", names(result))))
  expect_true(all(grepl("±", result[[grep("Estimate±SE", names(result))]])))
  expect_true(all(grepl("[0-9]", result[[grep("p-value", names(result))]])))
})

test_that("describe_emmeans handles asymptotic and regular CIs", {
  library(emmeans)
  model <- lm(Sepal.Length ~ Sepal.Width + Species, data = iris)
  emm_table <- emmeans(model, ~Species)
  # Should work for regular CIs
  expect_silent(rchiro::describe_emmeans(emm_table))
  # Should work for summary with infer=TRUE (asymptotic)
  expect_silent(rchiro::describe_emmeans(summary(emm_table, infer = c(TRUE, TRUE))))
})

# Edge case: missing df (should handle NA gracefully)
test_that("describe_contrasts and describe_emmeans handle missing df", {
  library(emmeans)
  model <- lm(Sepal.Length ~ Species, data = iris)
  contrast_table <- emmeans(model, pairwise ~ Species)$contrasts
  emm_table <- emmeans(model, ~Species)
  expect_silent(rchiro::describe_contrasts(contrast_table))
  expect_silent(rchiro::describe_emmeans(emm_table))
})
