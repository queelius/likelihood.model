

# Check if the print method is implemented for `class_name`
method_exists <- function(method, classname){
    
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
