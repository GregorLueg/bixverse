library(S7)

# Define your class
MyClass <- new_class(
  "MyClass",
  properties = list(
    name = class_character,
    value = class_numeric
  )
)

# Define S7 methods
method(format, MyClass) <- function(x, ...) {
  sprintf("MyClass object:\n  Name: %s\n  Value: %g", x@name, x@value)
}

method(str, MyClass) <- function(x, ...) {
  cat(format(x), "\n")
}

test = MyClass(name = "carl", value = 2)

class(test)

print(test)

test
