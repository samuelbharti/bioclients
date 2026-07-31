# biohttp's cache and circuit breaker are process-global singletons by design, so
# without this a test that trips the breaker or warms the cache silently changes
# the answer for whatever runs next.
withr::local_envvar(
  BIOHTTP_CACHE_DIR = withr::local_tempdir(.local_envir = teardown_env()),
  BIOHTTP_CACHE_SALT = "bioclients-test",
  BIOHTTP_CACHE_DISK = "",
  .local_envir = teardown_env()
)

biohttp::breaker_reset()
biohttp::cache_reset()

# Read a stored response body the way biohttp would hand it to a parser.
#
# simplifyVector = FALSE matters: biohttp's get_json() and post_json() parse with
# it off, so parsers walk nested lists rather than data frames. Reading a fixture
# any other way would test a shape the package never actually sees.
read_fixture <- function(name) {
  testthat::skip_if_not_installed("jsonlite")
  jsonlite::fromJSON(
    testthat::test_path("fixtures", name),
    simplifyVector = FALSE
  )
}

# A mocked JSON response, for the client halves.
mock_json <- function(body, status = 200L) {
  httr2::response(
    status_code = status,
    headers = list(`content-type` = "application/json"),
    body = charToRaw(body)
  )
}

# Every client test starts from a cold cache and a closed breaker.
reset_transport <- function() {
  biohttp::breaker_reset()
  biohttp::cache_reset()
}
