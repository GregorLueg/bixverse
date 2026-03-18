.onLoad <- function(...) {
  S7::methods_register()
  # S3 <> S7 weirdness
  registerS3method("print", "ScenicGrn", print.ScenicGrn)
}
