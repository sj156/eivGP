test_that('numerical CLI parsing preserves names, paths and literal values', {
  expect_identical(mixedgp_parse_cli_args(character()), list())
  expect_identical(mixedgp_parse_cli_args(c('--action=run', '--cores=32', '--data-root=/tmp/a=b c')),
    list(action = 'run', cores = '32', `data-root` = '/tmp/a=b c'))
  expect_identical(mixedgp_parse_cli_args('--methods=FULL,CUT')$methods, 'FULL,CUT')
})

test_that('numerical CLI parsing rejects ambiguous or missing values', {
  expect_error(mixedgp_parse_cli_args(c('--cores=4', '--cores=8')), 'Repeated option: cores')
  for (arg in c('cores=4', '--cores', '--cores=', ''))
    expect_error(mixedgp_parse_cli_args(arg), 'Use --name=value', fixed = TRUE)
})
