# Functions

rescale <- function(x) {
  (x-min(x))/(max(x)-min(x))
}
