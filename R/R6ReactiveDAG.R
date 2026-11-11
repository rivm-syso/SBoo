ReactiveDAG <- R6::R6Class(
  "ReactiveDAG",
  public = list(
    session = NULL,
    nodes = NULL,
    sources = NULL,
    params = NULL,
    debugR = NULL,
    
    initialize = function(params = list(), parent = NULL, debugR = FALSE) {
      # add parent (i.e. an SBcore) into params; it will be non-reactive, like other params
      if (!is.null(parent)) {
        params$parent <- parent
      }
      
      self$session <- shiny::MockShinySession$new()
      self$nodes   <- list()
      self$sources <- list()
      self$params  <- list2env(params, parent = emptyenv())
      self$debugR  <- debugR
    },
    
    add_source = function(name, value = NULL) {
      shiny::withReactiveDomain(self$session, {
        self$sources[[name]] <- shiny::reactiveVal(value)
      })
    },
    
    set_source = function(name, value) {
      if (!name %in% names(self$sources)) {
        self$add_source(name, value)
      } else {
        self$sources[[name]](value)
      }
    },
    
    add_node_reactive = function(node_name, node_fun_name = node_name, debug = self$debugR) {
      node_fun <- get(node_fun_name, mode = "function")
      
      wrapped <- private$make_smart_wrapper(
        fun          = node_fun,
        fun_name     = node_fun_name,
        fun_env      = parent.frame(),
        get_reactive = function(name) private$get_reactive(name),
        debug        = debug
      )
      
      shiny::withReactiveDomain(self$session, {
        self$nodes[[node_name]] <- shiny::reactive({
          wrapped()
        })
      })
    },
    
    get_value = function(name) {
      #browser()
      if (name %in% names(self$sources)) {
        shiny::isolate(self$sources[[name]]())
      } else if (name %in% names(self$nodes)) {
        shiny::isolate(self$nodes[[name]]())
      } else if (exists(name, envir = self$params, inherits = FALSE)) {
        get(name, envir = self$params, inherits = FALSE)
      } else {
        stop("Unknown dependency: ", name, call. = FALSE)
      }
    }
  ),
  
  active = list(
    knowns = function(value){
      if (missing(value)) {
        c(names(self$nodes), names(self$sources), names(self$params))
      } else { 
        stop("knowns is read-only")
      }
    },
    dataNames = function(value){
      if (missing(value)){
        names(self$sources)
      } else { 
        stop("dataNames is read-only")
      }
    }
  ),
  
  private = list(
    get_reactive = function(name) {
      lookup_name <- sub("^(to|all)\\.", "", name)
      
      if (lookup_name %in% names(self$sources)) {
        self$sources[[lookup_name]]()
      } else if (lookup_name %in% names(self$nodes)) {
        self$nodes[[lookup_name]]()
      } else if (exists(lookup_name, envir = self$params, inherits = FALSE)) {
        get(lookup_name, envir = self$params, inherits = FALSE)
      } else {
        stop("Unknown dependency: ", name, call. = FALSE)
      }
    },
    
    make_smart_wrapper = function(fun, fun_name, fun_env, get_reactive, debug = FALSE) {
      fun_formals <- names(formals(fun))
      arg_names <- fun_formals[fun_formals != "..."]
      
      function() {
        args <- setNames(lapply(arg_names, get_reactive), arg_names)
        
        if (debug) {
          message("Running node: ", fun_name)
        }
        
        do.call(fun, args)
      }
    }
  )
)