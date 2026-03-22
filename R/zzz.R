.onLoad <- function(...) {
  S7::methods_register()
  # S3 <> S7 weirdness
  registerS3method("print", "ScenicGrn", print.ScenicGrn)
  registerS3method("print", "CellQc", print.CellQc)
  registerS3method("print", "ScrubletRes", print.ScrubletRes)
  registerS3method("print", "BoostRes", print.BoostRes)
}
