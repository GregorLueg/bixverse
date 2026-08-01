.onLoad <- function(...) {
  S7::methods_register()

  # S3 <> S7 weirdness: manual registration of print methods
  classes <- c(
    "ScenicGrn",
    "CellQc",
    "ScrubletRes",
    "BoostRes",
    "ScDblFinderRes",
    "SingleCellNearestNeighbour",
    "KbetScores",
    "BatchSilhouetteScores",
    "BatchLisiScores",
    "Hotspot",
    "SingleCellFastClusters",
    "ADTCounts",
    "NmfResult",
    "StabilisedNmfResult",
    "LigandTargetInfluence",
    "ScTypeResults",
    "ScTypeCellResults",
    "BulkModuleResult",
    "ScPipeline",
    "ScStep",
    "SpCache",
    "SpatialSample",
    "SpMoransRes",
    "SpSparkxRes",
    "SpNhoodRes",
    "SpImageFeatures"
  )
  for (cls in classes) {
    registerS3method("print", cls, get(paste0("print.", cls)))
  }

  registerS3method("format", "BulkModuleResult", format.BulkModuleResult)
  registerS3method("dim", "BulkModuleResult", dim.BulkModuleResult)

  # same S3 <> S7 story for the plot methods on the spatial result classes
  for (cls in c("SpMoransRes", "SpSparkxRes", "SpNhoodRes")) {
    registerS3method("plot", cls, get(paste0("plot.", cls)))
  }
}
