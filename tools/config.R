# Note: Any variables prefixed with `.` are used for text
# replacement in the Makevars.in and Makevars.win.in

# check the packages MSRV first
source("tools/msrv.R")

# check DEBUG, DEV_BUILD and NOT_CRAN environment variables
env_debug <- Sys.getenv("DEBUG")
env_dev <- Sys.getenv("DEV_BUILD")
env_not_cran <- Sys.getenv("NOT_CRAN")

# check if the vendored zip file exists
vendor_exists <- file.exists("src/rust/vendor.tar.xz")

is_not_cran <- env_not_cran != ""
is_debug <- env_debug != ""

# DEV_BUILD keeps the cargo registry and the target dir warm between installs.
# It must not engage on a vendored build: those need CARGO_HOME to stay local
# so that vendor-config.toml is picked up.
is_dev <- env_dev != "" && !vendor_exists && !dir.exists("src/vendor")

if (is_debug) {
  # if we have DEBUG then we set not cran to true
  # CRAN is always release build
  is_not_cran <- TRUE
  message(
    "DEBUG requested but ignored - this package always builds release. ",
    "Set DEV_BUILD instead for a faster release build."
  )
}
is_debug <- FALSE

if (is_dev) {
  is_not_cran <- TRUE
  message(
    "DEV_BUILD set: keeping `src/rust/target` and using the shared cargo ",
    "registry. Still a release build, but without LTO - do not use this for ",
    "benchmarking or for a submission."
  )
}

if (!is_not_cran) {
  message("Building for CRAN.")
}

# we set cran flags only if NOT_CRAN is empty and if
# the vendored crates are present.
.cran_flags <- ifelse(
  !is_not_cran && vendor_exists,
  "-j 2 --offline",
  ""
)

# when DEBUG env var is present we use `--debug` build
.profile <- ifelse(is_debug, "", "--release")
.clean_targets <- ifelse(is_debug || is_dev, "", "$(TARGET_DIR)")

# used to replace @CARGO_HOME@. A CRAN build must not write outside the package,
# so cargo home points at a throwaway `src/.cargo`, which means a cold registry
# on every single install. Dev builds use the real one instead.
#
# The `rm -Rf` calls in Makevars deliberately reference CARGOTMP and not this,
# or a dev build would delete the user's cargo home.
.cargo_home <- ifelse(is_dev, "$(HOME)/.cargo", "$(CARGOTMP)")

# used to replace @DEV_EXPORTS@. Profile overrides come from the environment
# rather than a `[profile.release-dev]` because Cargo.toml is not owned by this
# template. opt-level stays at whatever the package declares so testing on real
# data remains feasible; what goes is the LTO link over the whole dependency
# graph, plus the codegen-units and strip settings that only pay off in a
# shipped build.
.dev_exports <- ifelse(
  is_dev,
  paste0(
    "CARGO_INCREMENTAL=1 ",
    "CARGO_PROFILE_RELEASE_LTO=false ",
    "CARGO_PROFILE_RELEASE_CODEGEN_UNITS=16 ",
    "CARGO_PROFILE_RELEASE_STRIP=none "
  ),
  ""
)

# We specify this target when building for webR
webr_target <- "wasm32-unknown-emscripten"

# here we check if the platform we are building for is webr
is_wasm <- identical(R.version$platform, webr_target)

# print to terminal to inform we are building for webr
if (is_wasm) {
  message("Building for WebR")
}

# we check if we are making a debug build or not
# if so, the LIBDIR environment variable becomes:
# LIBDIR = $(TARGET_DIR)/{wasm32-unknown-emscripten}/debug
# this will be used to fill out the LIBDIR env var for Makevars.in
target_libpath <- if (is_wasm) "wasm32-unknown-emscripten" else NULL
cfg <- if (is_debug) "debug" else "release"

# used to replace @LIBDIR@
.libdir <- paste(c(target_libpath, cfg), collapse = "/")

# use this to replace @TARGET@
# we specify the target _only_ on webR
# there may be use cases later where this can be adapted or expanded
.target <- ifelse(is_wasm, paste0("--target=", webr_target), "")

# add panic exports only for WASM builds
.panic_exports <- ifelse(
  is_wasm,
  "CARGO_PROFILE_DEV_PANIC=\"abort\" CARGO_PROFILE_RELEASE_PANIC=\"abort\" ",
  ""
)

# read in the Makevars.in file checking
is_windows <- .Platform[["OS.type"]] == "windows"

# if windows we replace in the Makevars.win.in
mv_fp <- ifelse(
  is_windows,
  "src/Makevars.win.in",
  "src/Makevars.in"
)

# set the output file
mv_ofp <- ifelse(
  is_windows,
  "src/Makevars.win",
  "src/Makevars"
)

# delete the existing Makevars{.win/.wasm}
if (file.exists(mv_ofp)) {
  message("Cleaning previous `", mv_ofp, "`.")
  invisible(file.remove(mv_ofp))
}

# read as a single string
mv_txt <- readLines(mv_fp)

# replace placeholder values
new_txt <- gsub("@CRAN_FLAGS@", .cran_flags, mv_txt) |>
  gsub("@PROFILE@", .profile, x = _) |>
  gsub("@CLEAN_TARGET@", .clean_targets, x = _) |>
  gsub("@LIBDIR@", .libdir, x = _) |>
  gsub("@TARGET@", .target, x = _) |>
  gsub("@PANIC_EXPORTS@", .panic_exports, x = _) |>
  gsub("@CARGO_HOME@", .cargo_home, x = _) |>
  gsub("@DEV_EXPORTS@", .dev_exports, x = _)

message("Writing `", mv_ofp, "`.")
con <- file(mv_ofp, open = "wb")
writeLines(new_txt, con, sep = "\n")
close(con)

message("`tools/config.R` has finished.")
