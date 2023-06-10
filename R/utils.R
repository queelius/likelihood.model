

# Check if the print method is implemented for `class_name`
method_exists <- function(generic_method, class_name){
    # Get all methods for the generic method
    methods <- methods(generic_method)

    # Create a logical vector, where each element is TRUE if the method is
    # implemented for the class, and FALSE otherwise
    impl_methods <- grepl(paste0(generic_method, ".", class_name), methods)

    # Check if any of the elements of the logical vector are TRUE
    any(impl_methods)
}
