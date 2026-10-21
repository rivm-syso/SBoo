#' Create a dynamic combiner function for reactive inputs
#'
#' @description
#' Constructs a function whose formal arguments are given by \code{inputs},
#' and whose body computes the sum of all those arguments. In the context of
#' \code{ReactiveDAG}, these formal arguments will be turned into active
#' bindings, so the resulting function acts as a reactive node that combines
#' multiple other reactives. This is needed for combining kaas of other processes but molecular deposition
#'
#' The generated function has the form:
#' \preformatted{
#' function(a, b, c, ...) {
#'   vals <- mget(ls(), inherits = FALSE)
#'   sum(unlist(vals))
#' }
#' }
#' where \code{a}, \code{b}, \code{c} are taken from \code{inputs}.
#'
#' @param inputs A character vector of names (e.g., \code{c("a", "b", "c")})
#'   corresponding to reactive sources, nodes, or parameters that will be
#'   combined.
#'
#' @return A function object with formal arguments matching \code{inputs},
#'   whose body computes the sum of those arguments.
#'
#' @examples
#' # Suppose you have a ReactiveDAG with nodes/sources "a", "b", "c"
#' dep_nodes <- c("a", "b", "c")
#' comb_fun <- make_combiner_fun(dep_nodes)
#'
#' # Register it as a node in your ReactiveDAG (assuming 'dag' is a ReactiveDAG)
#' # sum_dynamic <- comb_fun
#' # dag$add_node_reactive("sum_dynamic", "sum_dynamic")
#'
#' @export
make_combiner_fun <- function(inputs) {
  # inputs: character vector like c("a", "b", "c")
  
  body_expr <- quote({
    # formal args (a, b, c, ...) exist as active bindings created by ReactiveDAG
    vals <- mget(ls(), inherits = FALSE)  # get all formal args from this env
    sum(unlist(vals))
  })
  
  f <- function() {}
  formals(f) <- setNames(vector("list", length(inputs)), inputs)
  body(f) <- body_expr
  
  f
}