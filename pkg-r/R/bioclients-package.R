#' @keywords internal
#'
#' @details
#' Nothing is exported yet. The package is scaffolded and waiting on `biohttp`
#' to become installable, which is what Phase 1 of the build plan is gated on.
#'
#' The shape every service module will take, once they land: a request function
#' that calls one of `biohttp`'s entry points and returns its envelope, and a
#' parser function that takes an already-parsed body and returns a canonical
#' structure. The parser never touches the network, so a caller that already has
#' a response body can use it alone and the offline tests exercise it directly.
#'
#' @section On the empty Imports field:
#' `biohttp`, `jsonlite`, `rlang`, and `tibble` all belong in `Imports` and are
#' deliberately not there yet. `biohttp` is not installable from any repository
#' until it reaches r-universe, and declaring it now would fail `R CMD check`
#' for everyone including CI. Declaring the other three without it would leave
#' three imports that nothing calls, which is its own check note. All four go in
#' together with the first client.
"_PACKAGE"
