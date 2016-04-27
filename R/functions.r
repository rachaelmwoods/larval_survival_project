# Functions

rescale <- function(x, var=x) {
  (x-min(var))/(max(var)-min(var))
}

scaleback <- function(x, var) {
  (x * (max(var) - min(var))) + min(var)
}

