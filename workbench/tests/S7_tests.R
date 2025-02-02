library(S7)

# Define an S7 class
Person <- new_class("Person",
                    properties = list(
                      name = class_character,
                      age = class_numeric
                    )
)

# Define a custom print method
print_person <- function(x, ...) {
  cat(format(x), "\n")  # Use format() to ensure consistency
}

# Define a format method (used when object is printed directly)
format_person <- function(x, ...) {
  paste0("Person Object:\n",
         "Name: ", x@name, "\n",
         "Age: ", x@age)
}

# Register the methods
method(print, Person) <- print_person
method(format, Person) <- format_person  # Ensures auto-printing works

# Create an instance
p <- Person(name = "Alice", age = 30)

# Print the object
p  # Now this should trigger the custom format
