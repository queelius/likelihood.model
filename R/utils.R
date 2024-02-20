# Check if the print method is implemented for `class_name`
method_exists <- function(method, classname) {
  # Get all methods for the generic method
  the_methods <- methods(method)

  # Create a logical vector, where each element is TRUE if the method is
  # implemented for the class, and FALSE otherwise, and return
  # TRUE if any element is TRUE
  any(grepl(paste0(method, ".", classname), the_methods))
}

# Retrieve default parameters. If `par` is NULL or all NA, return `default`.
# Otherwise, return `par` with `default` values substituted for NA values.
# This may be useful in the context of `optim` where we want to use the
# default parameters if the user does not specify any. However, if the user
# specifies some parameters, we want to use those instead of the default
# parameters.
get_params <- function(par, default = NULL) {
  if (is.null(par) || all(is.na(par))) {
    return(default)
  }
  if (!is.null(default)) {
    par[is.na(par)] <- default[is.na(par)]
  }
  par
}

prepare_args <- function(primary_input, par, dist_function) {
  # Identify control arguments that should not be automatically assigned from 'par'
  control_args <- c("lower.tail", "log.p", "log")

  # Prepare the initial args list with the primary input
  args_list <- list(primary_input)

  if (!is.null(names(par))) {
    args_list <- c(args_list, as.list(par))
  } else {
    param_names <- names(formals(dist_function)[-1])
    param_names <- param_names[!param_names %in% control_args]
    unnamed_args_list <- setNames(as.list(par), param_names)
    args_list <- c(args_list, unnamed_args_list)
  }

  args_list
}
