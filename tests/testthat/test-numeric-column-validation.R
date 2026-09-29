test_that("numeric validation preserves values and invisible returns", {
  check <- IntegMultiReg:::.imr_check_numeric_columns
  inputs <- list(
    data.frame(id = 1:3, integer = 1:3, real = c(1, 2, 3)),
    data.frame(id = integer(), value = numeric()),
    data.frame(id = 1:3),
    data.frame(id = 1:3, value = I(1:3)),
    data.frame(id = 1:3, block = I(matrix(1:6, 3)))
  )
  for (input in inputs) {
    result <- withVisible(check(input, "input"))
    expect_identical(result$value, input)
    expect_false(result$visible)
  }
})

test_that("finite validation covers each selected column without changing error priority", {
  check <- IntegMultiReg:::.imr_check_numeric_columns
  for (value in list(NA_integer_, NA_real_, NaN, Inf, -Inf)) {
    for (column in c("first", "last")) {
      input <- data.frame(id = 1:3, first = 1:3, last = 1:3)
      input[[column]][2L] <- value
      expect_error(check(input, "input"), "All non-id values in `input` must be finite.",
                    fixed = TRUE)
    }
  }
  input <- data.frame(id = 1:3, invalid = c(Inf, 1, 2), text = letters[1:3])
  expect_error(check(input, "input"), "problem column: `text`", fixed = TRUE)
  expect_identical(withVisible(check(input, "input", character()))$value, input)
  expect_error(check(input, "input", "invalid"), "must be finite")
  input$invalid <- 1:3
  expect_identical(withVisible(check(input, "input", "invalid"))$value, input)
})

test_that("classed numeric columns retain matrix-coercion validation", {
  check <- IntegMultiReg:::.imr_check_numeric_columns
  input <- data.frame(id = 1:3, value = I(c(1, Inf, 3)))
  expect_error(check(input, "input"), "must be finite")
  input <- data.frame(id = 1:3, block = I(matrix(c(1:5, NA_real_), 3)))
  expect_error(check(input, "input"), "must be finite")
  input <- data.frame(id = 1:3, value = factor(letters[1:3]))
  expect_error(check(input, "input"), "must be numeric")
  input <- data.frame(id = 1:3, value = complex(real = 1:3, imaginary = 3:1))
  expect_error(check(input, "input"), "must be numeric")
})
