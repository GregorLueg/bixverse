library(S7)

# Define an S7 class
Person <- new_class("Person",
                    properties = list(
                      name = class_character,
                      age = class_numeric
                    )
)

# Define a custom print method
print_person <- function(x, ...) {  # Change 'object' to 'x'
  cat("Person Object:\n")
  cat("Name:", x@name, "\n")
  cat("Age:", x@age, "\n")
}

# Register the method for print()
method(print, Person) <- print_person

# Create an instance
p <- Person(name = "Alice", age = 30)

# Print the object
print(p)  # This will call the custom print method
p  # Also triggers print() in the console
