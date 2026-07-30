library(R6)
library(shiny)

ReactiveDAG <- R6Class(
  "ReactiveDAG",
  public = list(
    session = NULL,
    nodes = NULL,
    sources = NULL,
    params = NULL,
    
    initialize = function(params = list()) {
      self$session <- shiny::MockShinySession$new()
      self$nodes <- list()
      self$sources <- list()
      self$params <- new.env(parent = emptyenv())
    },
    
    add_source = function(name, value = NULL) {
      shiny::withReactiveDomain(self$session, {
        self$sources[[name]] <- reactiveVal(value)
      })
    },
    
    set_source = function(name, value) {
      self$sources[[name]](value)
    },
    
    add_node_reactive = function(node_fun_name) {
      node_fun <- get(node_fun_name)
      shiny::withReactiveDomain(self$session, {
        self$nodes[[name]] <- reactive({
          fmls <- formals(node_fun)
          arg_names <- names(fmls)
          
          args <- lapply(arg_names, function(nm) {
            if (nm %in% names(self$nodes)) {
              self$nodes[[nm]]()
            } else if (nm %in% names(self$sources)) {
              self$sources[[nm]]()
            } else if (nm %in% names(self$params)) {
              get(nm, envir = self$param_env, inherits = FALSE)
            } else {
              eval(fmls[[nm]])
            }
          })
          names(args) <- arg_names
          
          do.call(node_fun, args)
        })
      })
    },
    
    trigger = function(name) {
      shiny::withReactiveDomain(self$session, {
        isolate(self$nodes[[name]]())
      })
    }
  )
)

