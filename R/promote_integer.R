# Establish double arithmetic for ordinary integer inputs without changing the
# treatment of factors, classed objects, or other storage types. Unlike
# as.double(), storage.mode<- retains names, dimensions, and other attributes.
promote_integer <- function(x) {
  if (is.integer(x) && !is.object(x)) {
    storage.mode(x) <- "double"
  }
  x
}
