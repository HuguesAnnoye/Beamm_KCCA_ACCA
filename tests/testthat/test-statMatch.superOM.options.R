test_that("create default options", {
  expect_success(expect_type(statMatch.superOM.options(), "list"))
  expect_success(expect_s3_class(statMatch.superOM.options(), "superOM.options"))
})
