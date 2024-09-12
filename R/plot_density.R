calc_CI <- function(vector, alpha) {

  se <- sd(vector)/sqrt(length(vector))
  t <- qt((1 - alpha)/2 + 0.5, length(vector) - 1)
  CI = t*se
  CI

}

calc_density_over_UMAP <- function(obj) {

  x <- obj@reductions$umap.cca@cell.embeddings[, 1]
  y <- obj@reductions$umap.cca@cell.embeddings[, 2]
  data <- FetchData(obj, assay = "RNA", layer = "data", vars = "CD276")[, "CD276"]

  z <- calc_density(x, y, data)
  names(z) <- colnames(obj)
  return(z)

}

calc_density <- function(x, y, data) {

  w <- data / sum(data) * length(data)
  h <- c(
    ks::hpi(x),
    ks::hpi(y)
  )
  h <- h * 0.1
  dens <- ks::kde(x = matrix(data = c(x, y), ncol = 2), w = w / sum(w) * length(w), bgridsize = rep(1000, 2))
  ix <- findInterval(x, dens$eval.points[[1]])
  iy <- findInterval(y, dens$eval.points[[2]])
  ii <- cbind(ix, iy)
  z <- dens$estimate[ii]
  return(z)

}

plot_umap <- function(obj, feature, feature_title) {

  x <- obj@reductions$umap.cca@cell.embeddings[, 1]
  y <- obj@reductions$umap.cca@cell.embeddings[, 2]
  z <- FetchData(obj, assay = "RNA", layer = "data", vars = feature)[, feature]
  obj$z <- calc_density(x, y, z)
  umap_plot <- FeaturePlot(obj,
                           features = c("z"),
                           raster = FALSE,
                           pt.size = 0.0005,
                           label = FALSE,
                           label.color = "black"
                           ) +
    ggtitle(label = features_title) +
    theme(plot.title = element_text(vjust = - 3.5),
          plot.margin = unit(c(-5, -5, -5, -5), units = "mm")
    ) +
    NoLegend() +
    NoAxes()

  return(umap_plot)

}
