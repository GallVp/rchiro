# Additional tests for describe_contrasts and describe_emmeans with lmer models

test_that("describe_contrasts works with lmer model", {
  library(lme4)
  library(emmeans)
  # Use cbpp dataset: size ~ period + (1 | herd)
  data(cbpp, package = "lme4")
  model <- lmer(size ~ period + (1 | herd), data = cbpp)
  contrast_table <- emmeans(model, pairwise ~ period)$contrasts
  result <- rchiro::describe_contrasts(contrast_table)
  expect_true(is.data.frame(result))
  expect_true(any(grepl("Difference±SE", names(result))))
  expect_true(any(grepl("p-value", names(result))))
})

test_that("describe_emmeans works with lmer model", {
  library(lme4)
  library(emmeans)
  # Use cbpp dataset: size ~ period + (1 | herd)
  data(cbpp, package = "lme4")
  model <- lmer(size ~ period + (1 | herd), data = cbpp)
  emm_table <- emmeans(model, ~period)
  result <- rchiro::describe_emmeans(emm_table)
  expect_true(is.data.frame(result))
  expect_true(any(grepl("Estimate±SE", names(result))))
  expect_true(any(grepl("p-value", names(result))))
})
