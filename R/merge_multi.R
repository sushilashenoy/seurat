#' Merge multiple Seurat objects
#'
#' @param objects list of Seurat objects, or a single Seurat object if group.by will be used
#' @param ... Extra parameters to MergeSeurat
#'
#' @return Returns Seurat object data with data combined from list of objects
#'
#' @seealso \code{MergeSeurat}
#'
#' @export
#'
#'
MergeMulti <- function(
  objects,
  ...
) {

  if ( missing(x = objects) || !is.list(objects) ) {
    stop("objects must be a list of seurat objects.")
  }
   
  if ( is.list(objects) && length(objects) > 1 ) {
    cat("Merging objects\n", file = stderr())
    combined.object <- Reduce(function (object1, object2) {
      MergeSeurat(
      object1 = object1,
      object2 = object2,
      ...) }, objects)
    # to improve, to pull the same normalization and scale params as previously used
    combined.object <- ScaleData(object = combined.object)
    combined.object@scale.data[is.na(x = combined.object@scale.data)] <- 0
    combined.object@var.genes <- genes.use
   
    return(object)
  }
}

