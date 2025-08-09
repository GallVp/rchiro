# Tests for confint.rlmerMod

test_that("confint.rlmerMod returns correct confidence intervals for all parameters", {
  set.seed(123)
  library(robustlmm)
  library(lme4)
  data(sleepstudy, package = "lme4")
  mod <- robustlmm::rlmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  ci <- confint.rlmerMod(mod)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), length(lme4::fixef(mod)))
  expect_equal(colnames(ci), scales::label_percent(accuracy = 0.1)(c(0.025, 0.975)))
})

test_that("confint.rlmerMod returns correct confidence interval for a single parameter", {
  set.seed(123)
  library(robustlmm)
  library(lme4)
  data(sleepstudy, package = "lme4")
  mod <- robustlmm::rlmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  ci <- confint.rlmerMod(mod, parm = "Days")
  expect_true(is.numeric(ci))
  expect_equal(length(ci), 2)
  expect_equal(names(ci), scales::label_percent(accuracy = 0.1)(c(0.025, 0.975)))
})
