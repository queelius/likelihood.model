

# Check if the print method is implemented for `class_name`
method_exists <- function(method, classname){
    
    # Get all methods for the generic method
    the_methods <- methods(method)

    # Create a logical vector, where each element is TRUE if the method is
    # implemented for the class, and FALSE otherwise, and return
    # TRUE if any element is TRUE
    any(grepl(paste0(method, ".", classname), the_methods))
}
