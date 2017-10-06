#' Perform Canonical Correlation Analysis
#'
#' Runs a canonical correlation analysis using a diagonal implementation of CCA.
#' For details about stored CCA calculation parameters, see
#' \code{PrintCCAParams}.
#'
#' @param objects list of Seurat objects, or a single Seurat object if group.by will be used
#' @param group.by Factor to group by (column vector stored in object@@meta.data)
#' @param num.cc Number of canonical vectors to calculate
#' @param genes.use Set of genes to use in CCA. Default is object@@var.genes. If
#' two objects are given, the default is the union of both variable gene sets
#' that are also present in both objects.
#' @param scale.data Use the scaled data from the object
#' @param rescale.groups Rescale each set of cells independently
#' @param rgcca.tol tol parameter for RGCCA::rgcca()
#' @param ... Extra parameters to MergeSeurat
#'
#' @return Returns Seurat object with the CCA stored in the @@dr$cca slot. If
#' one object is passed, the same object is returned. If two or more are passed, a
#' combined object is returned.
#'
#' @seealso \code{MergeSeurat} \code{RunCCA}
#'
#' @importFrom RGCCA rgcca
#' @export
#'
#' @examples
#' pbmc_small
#' # As GCCA works with more than two datasets, we will split our test object into two just for this example
#' pbmc1 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[1:30])
#' pbmc2 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[31:60])
#' pbmc3 <- SubsetData(pbmc_small,cells.use = pbmc_small@cell.names[61:90])
#' pbmc1@meta.data$group <- "group1"
#' pbmc2@meta.data$group <- "group2"
#' pbmc3@meta.data$group <- "group3"
#' pbmc_cca <- RunGCCA(list(pbmc1,pbmc2,pbmc3))
#' # Print results
#' PrintDim(pbmc_cca,reduction.type = 'gcca')
#'
RunGCCA <- function(
  objects,
  group.by,
  num.cc = 20,
  genes.use,
  scale.data = TRUE,
  rescale.groups = FALSE,
  rgcca.tol=1e-3,
  ...
) {

  if ( missing(x = objects) ) {
    stop("objects must be a list of seurat objects or a single seurat object (with group.by defined).")
  }
  if (! missing(x = group.by) && ! missing(x = objects) && is.list(objects) && length(objects) > 1 ) {
    warning("Multiple objects and group.by set. Continuing with objects defining the groups")
  }
  if ( is.list(objects) && length(objects) > 1) {
    if (missing(x = genes.use) ) {
      genes.use <- Reduce(union, lapply(objects, function (x) x@var.genes))
      if (length(x = genes.use) == 0) {
        stop("No variable genes present. Run MeanVarPlot and retry")
      }
    }
    if (scale.data) {
      possible.genes <- Reduce(intersect, lapply(objects, function (x) rownames(x@scale.data)))
      genes.use <- genes.use[genes.use %in% possible.genes]
      data.use <- lapply(objects, function (x) x@scale.data[genes.use, ])
    } else {
      possible.genes <- Reduce(intersect, lapply(objects, function (x) rownames(x@data)))
      genes.use <- genes.use[genes.use %in% possible.genes]
      data.use <- lapply(objects, function (x) x@data[genes.use, ])
    }
    if (length(x = genes.use) == 0) {
      stop("0 valid genes in genes.use")
    }
  } else {
    if (! missing(x = group.by)) {
      if (! group.by %in% colnames(x = object@meta.data)) {
        stop("invalid group.by parameter")
      }
    }
    if (missing(x = genes.use)) {
      genes.use <- object@var.genes
      if (length(x = genes.use) == 0) {
        stop("No variable genes present. Run MeanVarPlot and retry")
      }
    }
    if (! missing(x = group.by)) {
      if ( is.list(objects) ) {
        object <- objects[[1]]
      } else {
        object <- objects
      }

      groups <- unique(object@meta.data[, group.by])

      object.current.ids <- object@ident
      object <- SetAllIdent(object = object, id = group.by)
      cells <- lapply(groups, function (x) CheckGroup(object, group = x, group.id=x))
      object <- SetIdent(
        object = object,
        cells.use = object@cell.names,
        ident.use = object.current.ids
      )
    }
    if (scale.data) {
      if (rescale.groups) {
        data.use <- lapply(1:length(groups), function (i) ScaleData(
          object = object,
          data.use = object@data[genes.use, cells[[i]]]
        )@scale.data)
        
      } else {
        data.use <- lapply(cells, function (x) object@scale.data[genes.use, x])
      }
    } else {
        data.use <- lapply(cells, function (x) object@data[genes.use, x])
    }
  }

  ngroups <- length(data.use)

  for ( i in 1:ngroups ) {
    genes.use <- CheckGenes(data.use = data.use[[i]], genes.use = genes.use)
  }

  data.use <- lapply(data.use, function (x) x[genes.use, ])

  cat("Running GCCA\n", file = stderr())

  gcca.results <- GenCanCor( mats = data.use, standardize = TRUE, k = num.cc, tol=rgcca.tol)
  gcca.data <- gcca.results$components
  colnames(x = gcca.data) <- paste0("GCC", 1:num.cc)
  if ( is.list(objects) && length(objects) > 1 ) {
    cat("Merging objects\n", file = stderr())
    combined.object <- Reduce(function (object1, object2) {
      MergeSeurat(
      object1 = object1,
      object2 = object2,
      do.scale = FALSE,
      do.center = FALSE,
      ...) }, objects)
    # to improve, to pull the same normalization and scale params as previously used
    combined.object <- ScaleData(object = combined.object)
    combined.object@scale.data[is.na(x = combined.object@scale.data)] <- 0
    combined.object@var.genes <- genes.use
    rownames(cca.data) <- colnames(combined.object@data)
    combined.object <- SetDimReduction(
      object = combined.object,
      reduction.type = "gcca",
      slot = "cell.embeddings",
      new.data = gcca.data
    )
    combined.object <- SetDimReduction(
      object = combined.object,
      reduction.type = "gcca",
      slot = "key",
      new.data = "GCC"
    )
    combined.object <- ProjectDim(
      object = combined.object,
      reduction.type = "gcca",
      do.print = FALSE
    )
    combined.object <- SetDimReduction(
      object = combined.object,
      reduction.type = "gcca",
      slot = "gene.loadings",
      new.data = GetGeneLoadings(
        object = combined.object,
        reduction.type = "gcca",
        use.full = TRUE,
        genes.use = genes.use
      )
    )
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunGCCA"))]
    combined.object <- SetCalcParams(
      object = combined.object,
      calculation = "RunGCCA",
      ... = parameters.to.store
    )
    return(combined.object)
  } else {
    object <- SetDimReduction(
      object = object,
      reduction.type = "gcca",
      slot = "cell.embeddings",
      new.data = gcca.data
    )
    object <- SetDimReduction(
      object = object,
      reduction.type = "gcca",
      slot = "key",
      new.data = "GCC"
    )

    object <- ProjectDim(
      object = object,
      reduction.type = "gcca",
      do.print = FALSE
    )
    object@scale.data[is.na(x = object@scale.data)] <- 0
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunGCCA"))]
    object <- SetCalcParams(
      object = object,
      calculation = "RunGCCA",
      ... = parameters.to.store
    )
    return(object)
  }
}




# Run the diagonal canonical correlation procedure
#
# @param mat1         First matrix
# @param mat2         Second matrix
# @param standardize  Standardize matrices - scales columns to have unit
#                     variance and mean 0
# @param k            Number of canonical correlation vectors (CCs) to calculate
#
# @return             Returns the canonical correlation vectors - corresponding
#                     to the left and right singular vectors after SVD - as well
#                     as the singular values.
#
GenCanCor <- function(mats, standardize = TRUE, k = 20, tol=1e-3) {
  set.seed(seed = 42)
  if (standardize) {
    mats <- lapply(mats, function (x) Standardize(mat = x, display_progress = FALSE))
  }
  nmat <- length(mats)
  mmat <- do.call(cbind, mats)

  C <- matrix(0, nmat+1, nmat+1)
  C[nmat+1, 1:nmat] <- C[1:nmat, nmat+1] <- 1

  mats[[nmat+1]] <- mmat

  gcca <- rgcca(A=mats, C=C, tau=rep(0, nmat+1),
    ncomp=c(rep(1, nmat), k), tol=tol, scheme='factorial')
  gcca$components <- gcca$a[[nmat+1]]

  return(gcca)
}
