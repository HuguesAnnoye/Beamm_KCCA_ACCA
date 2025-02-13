set.seed(1)
N <- 200
df.rec <- data.frame(
  CV_1 = runif(N),
  CV_2 = runif(N),
  NCV_rec_1 = runif(N),
  NCV_rec_2 = runif(N)
)
df.don <- data.frame(
  CV_1 = runif(N),
  CV_2 = runif(N),
  NCV_don_1 = runif(N),
  NCV_don_2 = runif(N)
)

test_that("check errors only continuous variables",{
  # Options
  opts <- "test"
  expect_error(statMatch(df.rec, df.don, opts = opts, method = "superOM"), "Argument 'opts' should be an object of class 'superOM.options'", fixed = TRUE)
  # Zero constraints
  expect_error(statMatch(df.rec, df.don, zero.constraints = c("NCV_don_1", "NCV_don_2"),method = "superOM"), "Problem in 'zero.constraints', no zeros in variable(s): NCV_don_1, NCV_don_2", fixed = TRUE)
  expect_error(statMatch(df.rec, df.don, zero.constraints = c("wrong_name1", "wrong_name2"),method = "superOM"), "Variable(s) in 'zero.constraints' not in 'names.NCV': wrong_name1, wrong_name2", fixed = TRUE)
})


