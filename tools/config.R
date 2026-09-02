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
.clean_targets <- ifelse(is_debug || is_dev, "", "\"$(TARGET_DIR)\"")

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

# used to replace @TARGET_DIR@. Windows only, and it is about MAX_PATH rather
# than taste. R CMD INSTALL builds under
# `AppData/Local/Temp/RtmpXXXXXX/R.INSTALLXXXXXXXXXX/bixverse/src/`, 84
# characters before the target dir starts. `hdf5-metno-src` is the one
# dependency that shells out to CMake, and CMake's TryCompile scratch adds a
# further 139, which puts the object path at 265 against a 260 limit. The
# symptom is not a path error, it is `gcc.exe` being declared "not able to
# compile a simple test program", because the .obj silently never lands.
#
# Building outside the package tree is what keeps it short. Kept out of
# `Makevars.in`: no unix path is anywhere near the limit.
.target_dir <- if (is_windows) {
  home <- Sys.getenv("USERPROFILE", unset = Sys.getenv("HOME"))
  d <- file.path(home, ".bixverse-cargo")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  normalizePath(d, winslash = "/", mustWork = TRUE)
} else {
  "./rust/target"
}

# HDF5 provider. `bixverse-rs` 0.5.0 made the from-source HDF5 build opt-in
# (`hdf5-static`, which this package turns on by default). Building it is the
# only option that needs no system library, and it is the right one everywhere
# except a cross-ABI Windows build: cargo running from an msvc host against a
# `-pc-windows-gnu` target makes `hdf5-metno-src` name its output the msvc way,
# and the link then fails looking for `libhdf5`. That is what the r-universe
# Windows jobs hit, and only them.
#
# There we look for an external libhdf5 through pkg-config instead.
# `hdf5-metno-sys` skips its runtime version check precisely when pkg-config
# wins on Windows, which is what makes a static-only Rtools HDF5 usable, and
# pkg-config also hands back the transitive link flags so nothing has to be
# guessed.
#
# The gate matters. A gnu host, which is what our own CI installs, builds HDF5
# from source correctly, and putting Rtools' lib dir on the link line there
# instead drags a second mingw runtime into it: the runner carries its own
# mingw, whose `crt2.o` then resolves `libmsvcrt.a` out of Rtools and loses
# `__p___initenv`.
.cargo_features <- ""
.hdf5_exports <- ""
.hdf5_libs <- ""
.hdf5_rustflags <- ""

# `rustc -vV` reports the toolchain's own triple. The target is always a
# `-pc-windows-gnu*` one here, since extendr links against the gnu ABI, so an
# msvc host is exactly the cross-ABI case.
rustc_host <- if (is_windows) {
  out <- tryCatch(
    system2("rustc", "-vV", stdout = TRUE, stderr = FALSE),
    error = function(e) character(0),
    warning = function(w) character(0)
  )
  host <- grep("^host: ", out, value = TRUE)
  if (length(host)) sub("^host: ", "", host[1]) else ""
} else {
  ""
}

is_cross_abi <- grepl("windows-msvc", rustc_host, fixed = TRUE)

if (is_windows && !is_cross_abi) {
  message(
    "Rust host `",
    rustc_host,
    "` matches the target ABI. Building HDF5 from source."
  )
}

if (is_windows && is_cross_abi) {
  pkg_config <- Sys.which("pkg-config")

  candidates <- character(0)
  if (nzchar(pkg_config)) {
    rtools_homes <- Sys.getenv(c(
      "RTOOLS45_AARCH64_HOME",
      "RTOOLS45_HOME",
      "RTOOLS44_HOME",
      "RTOOLS43_HOME",
      "RTOOLS42_HOME"
    ))
    rtools_homes <- unique(rtools_homes[nzchar(rtools_homes)])
    prefixes <- c(
      "x86_64-w64-mingw32.static.posix",
      "aarch64-w64-mingw32.static.posix",
      "clang-aarch64",
      "ucrt64",
      "mingw64"
    )
    candidates <- file.path(
      rep(rtools_homes, each = length(prefixes)),
      prefixes,
      "lib",
      "pkgconfig"
    )
    candidates <- candidates[dir.exists(candidates)]
    # whatever the caller already set stays first in line
    candidates <- c(Sys.getenv("PKG_CONFIG_PATH"), candidates)
    candidates <- unique(candidates[nzchar(candidates)])
  }

  # `system2(env = )` is documented as unsupported on Windows, which is the
  # only platform this branch runs on, so the variable is set and restored
  # around the probe instead.
  old_pkg_config_path <- Sys.getenv("PKG_CONFIG_PATH", unset = NA)
  on.exit(
    if (is.na(old_pkg_config_path)) {
      Sys.unsetenv("PKG_CONFIG_PATH")
    } else {
      Sys.setenv(PKG_CONFIG_PATH = old_pkg_config_path)
    },
    add = TRUE
  )

  for (path in candidates) {
    # forward slashes, or the Makevars recipe and pkg-config disagree about
    # what a backslash means
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    Sys.setenv(PKG_CONFIG_PATH = path)
    libs <- suppressWarnings(system2(
      pkg_config,
      c("--libs", "--static", "hdf5"),
      stdout = TRUE,
      stderr = FALSE
    ))
    if (!is.null(attr(libs, "status")) || !length(libs)) {
      next
    }
    .cargo_features <- "--no-default-features"
    # the recipe runs under sh, not cmd, whatever the host is
    .hdf5_exports <- paste0(
      "PKG_CONFIG_PATH=",
      shQuote(path, type = "sh"),
      " "
    )
    flags <- trimws(paste(libs, collapse = " "))
    # R links the final DLL with PKG_LIBS, but cargo links the `document`
    # binary itself, so the same flags have to go both ways. `hdf5-metno-sys`
    # emits a bare `-lhdf5` and drops pkg-config's Libs.private, which is why
    # zlib and szip have to be handed over explicitly. rustc takes `-L`/`-l` in
    # the attached form pkg-config already prints.
    .hdf5_libs <- paste0(flags, " ")
    .hdf5_rustflags <- paste0(" ", flags)
    message(
      "Found an external HDF5 via `",
      path,
      "`. Skipping the source build."
    )
    break
  }

  if (!nzchar(.cargo_features)) {
    message("No external HDF5 found. Building it from source.")
  }
}

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
  gsub("@DEV_EXPORTS@", .dev_exports, x = _) |>
  gsub("@TARGET_DIR@", .target_dir, x = _) |>
  # fixed = TRUE: these carry Windows paths, and a backslash in a gsub
  # replacement is an escape rather than a literal.
  gsub("@CARGO_FEATURES@", .cargo_features, x = _, fixed = TRUE) |>
  gsub("@HDF5_EXPORTS@", .hdf5_exports, x = _, fixed = TRUE) |>
  gsub("@HDF5_LIBS@", .hdf5_libs, x = _, fixed = TRUE) |>
  gsub("@HDF5_RUSTFLAGS@", .hdf5_rustflags, x = _, fixed = TRUE)

message("Writing `", mv_ofp, "`.")
con <- file(mv_ofp, open = "wb")
writeLines(new_txt, con, sep = "\n")
close(con)

message("`tools/config.R` has finished.")
