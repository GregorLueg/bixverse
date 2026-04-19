.onLoad <- function(...) {
  S7::methods_register()

  # S3 <> S7 weirdness
  # manual registering here
  registerS3method("print", "ScenicGrn", print.ScenicGrn)
  registerS3method("print", "CellQc", print.CellQc)
  registerS3method("print", "ScrubletRes", print.ScrubletRes)
  registerS3method("print", "BoostRes", print.BoostRes)
  registerS3method("print", "ScDblFinderRes", print.ScDblFinderRes)
  registerS3method(
    "print",
    "SingleCellNearestNeighbour",
    print.SingleCellNearestNeighbour
  )
  registerS3method(
    "print",
    "KbetScores",
    print.KbetScores
  )
  registerS3method(
    "print",
    "BatchSilhouetteScores",
    print.BatchSilhouetteScores
  )
  registerS3method(
    "print",
    "BatchLisiScores",
    print.BatchLisiScores
  )
  registerS3method(
    "print",
    "Hotspot",
    print.Hotspot
  )
}
