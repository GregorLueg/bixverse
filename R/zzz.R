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
    "ScTypeResults"
  )
  for (cls in classes) {
    registerS3method("print", cls, get(paste0("print.", cls)))
  }
}
