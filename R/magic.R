
`%^%` <- function (x, n) {
  if (n == 1 ) return(x)
  y <- x
  for ( i in 2:n ) {
    y <- x %*% y
  }
  return(y)
}

#' Impute expression based on neighboring cell using method based on
#' 
#' MAGIC: A diffusion-based imputation method reveals gene-gene interactions
#' in single-cell RNA-sequencing data
#' 
#' David van Dijk, Juozas Nainys, Roshan Sharma, Pooja Kathail, Ambrose J Carr,
#' Kevin R Moon, Linas Mazutis, Guy Wolf, Smita Krishnaswamy, Dana Pe'er
#' doi: https://doi.org/10.1101/111591
#' 
#' 
#' 
#' @export

RunMagic <- function(object, t=6) {
  if (.hasSlot(object = object, name = "snn")) {
    if (length(x = object@snn) > 1) {
      snn.built <- TRUE
    }
  }

  if ( !snn.built ) {
  	stop('Please run FindClusters with save.SNN=TRUE before running this function.')
  }

  if ( t != abs(round(t)) ) {
  	stop('t must be a positive integer.')
  }

  affinities <- object@snn + t(object@snn)

  message('Row normalizing affinity matrix..')
  markov <- sweep(affinities, 1, apply(affinities, 1, sum), '/')

  message('Calculating diffusion.. (t=', t, ')')
  diffusion <- markov %^% t

  message('Calculating imputed values')
  imp.data <- t(diffusion %*% t(dds@data))

  message('Calculating scaling factors')
  scaling.factor <- apply(dds@data, 1, quantile, 0.99)/apply(imp.data, 1, max)

  message('Saving scaled imputed data')
  scaled.imp.data <- sweep(imp.data, 1, scaling.factor, '*')
  object@imputed <- scaled.imp.data

  return ( object )
}
