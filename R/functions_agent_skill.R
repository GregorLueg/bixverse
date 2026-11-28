# agent skill -----------------------------------------------------------------

## helpers ---------------------------------------------------------------------

#' The skills directory an agent discovers skills in
#'
#' @param agent String. One of `c("claude", "codex")`.
#' @param scope String. One of `c("user", "project")`.
#'
#' @returns String. The directory, without the `bixverse` subdirectory. Project
#' paths are relative to the working directory.
#'
#' @keywords internal
.agent_skill_dir <- function(agent, scope) {
  # Codex scans the open standard's `.agents/skills`, not `~/.codex/skills`
  base <- switch(agent, claude = ".claude/skills", codex = ".agents/skills")
  if (scope == "user") file.path("~", base) else base
}

#' Whether paths are symbolic links
#'
#' @param paths Character. The paths.
#'
#' @returns Logical, `FALSE` for paths that do not exist.
#'
#' @keywords internal
.is_symlink <- function(paths) {
  # Sys.readlink() gives NA for a missing path and "" for a non-link
  links <- Sys.readlink(paths)
  !is.na(links) & nzchar(links)
}

#' Read the package version an installed skill was copied from
#'
#' @param target String. The installed skill directory.
#'
#' @returns String, or `NA` for an install without a version stamp.
#'
#' @keywords internal
.agent_skill_version <- function(target) {
  stamp <- file.path(target, ".bixverse-version")
  if (!file.exists(stamp)) {
    return(NA_character_)
  }
  trimws(readLines(stamp, n = 1L, warn = FALSE))
}

#' Warn about installed skills that lag the package
#'
#' @description
#' Called on attach. A copied skill does not follow a package upgrade, so an
#' agent would keep reading the API of the old version. Never errors.
#'
#' @returns `NULL`, invisibly.
#'
#' @keywords internal
.warn_stale_agent_skills <- function() {
  status <- tryCatch(agent_skill_status(), error = function(e) NULL)
  if (is.null(status)) {
    return(invisible(NULL))
  }
  for (i in which(!status$current)) {
    packageStartupMessage(sprintf(
      paste(
        "The bixverse agent skill at %s is from version %s, the package is %s.",
        "Run `bixverse::install_agent_skill(overwrite = TRUE)` to update it."
      ),
      status$path[i],
      if (is.na(status$version[i])) "unknown" else status$version[i],
      utils::packageVersion("bixverse")
    ))
  }
  invisible(NULL)
}

## installation ----------------------------------------------------------------

#' Install the bixverse agent skill
#'
#' @description
#' bixverse ships a skill that teaches coding agents how to *use* the package:
#' the workflows, the `params_*()` convention, the ordering constraints and the
#' traps. This copies it out of the installed package to where your agent looks
#' for it.
#'
#' The skill is versioned with the package. The install is stamped with the
#' package version, and `library(bixverse)` warns when an installed skill is
#' older than the package; re-run this with `overwrite = TRUE` after upgrading.
#' From a shell, see [bixverse::skill_cli()].
#'
#' @details
#' `SKILL.md` follows the Agent Skills open standard, so the same files serve
#' every agent; only the directory differs. `"claude"` installs to
#' `.claude/skills`, `"codex"` to `.agents/skills`, which Codex and other agents
#' scan. A project install ties the skill to the project, so two projects on
#' different bixverse versions each see their own API. Run it from the project
#' root.
#'
#' @param dest String or `NULL`. Skills directory to install into; the skill
#' lands in a `bixverse` subdirectory. Overrides `agent` and `scope`.
#' @param agent String. One of `c("claude", "codex")`.
#' @param scope String. One of `c("user", "project")`: under the home
#' directory, or under the working directory.
#' @param overwrite Boolean. Replace an existing installation. Defaults to
#' `FALSE`, in which case an existing `bixverse` directory causes an error
#' rather than a silent clobber.
#'
#' @returns The path the skill was written to, invisibly.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Claude Code, into ~/.claude/skills
#' install_agent_skill()
#'
#' # Codex and other open standard agents, into ./.agents/skills
#' install_agent_skill(agent = "codex", scope = "project")
#'
#' # refresh after a package upgrade
#' install_agent_skill(overwrite = TRUE)
#' }
install_agent_skill <- function(
  dest = NULL,
  agent = c("claude", "codex"),
  scope = c("user", "project"),
  overwrite = FALSE
) {
  agent <- match.arg(agent)
  scope <- match.arg(scope)

  # checks
  checkmate::qassert(dest, c("S1", "0"))
  checkmate::assertChoice(agent, c("claude", "codex"))
  checkmate::assertChoice(scope, c("user", "project"))
  checkmate::qassert(overwrite, "B1")

  source_dir <- system.file("skills", "bixverse", package = "bixverse")
  if (source_dir == "") {
    stop(paste(
      "The agent skill was not found in the installed package. This usually",
      "means bixverse was installed from a source tree predating the skill."
    ))
  }

  dest <- path.expand(
    if (is.null(dest)) .agent_skill_dir(agent, scope) else dest
  )
  target <- file.path(dest, "bixverse")

  if (dir.exists(target) || .is_symlink(target)) {
    if (!overwrite) {
      stop(sprintf(
        "A bixverse skill already exists at %s. Re-run with `overwrite = TRUE`.",
        target
      ))
    }
    unlink(target, recursive = TRUE, force = TRUE)
  }

  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  # copies source_dir itself (not just its contents), giving dest/bixverse
  ok <- file.copy(source_dir, dest, recursive = TRUE)
  if (!ok || !file.exists(file.path(target, "SKILL.md"))) {
    stop(sprintf("Failed to install the agent skill to %s", target))
  }
  writeLines(
    as.character(utils::packageVersion("bixverse")),
    file.path(target, ".bixverse-version")
  )

  message(sprintf("bixverse agent skill installed to %s", target))
  message("Start a new agent session to pick it up.")

  invisible(target)
}

#' Installed bixverse agent skills and their versions
#'
#' @description
#' Looks in the user and project locations of every supported agent. Symlinked
#' installs track their source, so they are left out.
#'
#' @returns A data.table with `path`, `version` (the stamped package version,
#' `NA` for an install without a stamp) and `current` (whether it matches the
#' installed package).
#'
#' @export
#'
#' @examples
#' agent_skill_status()
agent_skill_status <- function() {
  dirs <- c(
    .agent_skill_dir("claude", "user"),
    .agent_skill_dir("claude", "project"),
    .agent_skill_dir("codex", "user"),
    .agent_skill_dir("codex", "project")
  )
  targets <- file.path(path.expand(dirs), "bixverse")
  targets <- targets[dir.exists(targets) & !.is_symlink(targets)]

  versions <- vapply(targets, .agent_skill_version, character(1L))
  data.table::data.table(
    path = targets,
    version = unname(versions),
    current = !is.na(versions) &
      versions == as.character(utils::packageVersion("bixverse"))
  )
}

## command line ----------------------------------------------------------------

#' Manage the agent skill from a shell
#'
#' @description
#' The command line of the agent skill, the R twin of the Python
#' `bixverse skill` command:
#'
#' ```
#' Rscript -e 'bixverse::skill_cli()' skill install [--agent claude|codex]
#'   [--project] [--dest DIR] [--force]
#' Rscript -e 'bixverse::skill_cli()' skill status
#' ```
#'
#' Run it with the R library the project uses, so the skill matches the
#' bixverse version the project sees.
#'
#' @param args Character. The command line arguments. Defaults to the trailing
#' arguments of the R session.
#'
#' @returns The exit code, invisibly: `1` from `status` if any install lags the
#' package, else `0`. Errors on bad arguments.
#'
#' @export
#'
#' @examples
#' skill_cli(c("skill", "status"))
skill_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  checkmate::qassert(args, "S*")

  usage <- paste(
    "Usage:",
    "  Rscript -e 'bixverse::skill_cli()' skill install [--agent claude|codex]",
    "    [--project] [--dest DIR] [--force]",
    "  Rscript -e 'bixverse::skill_cli()' skill status",
    sep = "\n"
  )

  if (length(args) == 0L || any(args %in% c("-h", "--help"))) {
    cat(usage, "\n")
    return(invisible(0L))
  }
  if (identical(args, "--version")) {
    cat(as.character(utils::packageVersion("bixverse")), "\n")
    return(invisible(0L))
  }
  if (length(args) < 2L || args[1L] != "skill") {
    stop(usage, call. = FALSE)
  }

  action <- args[2L]
  flags <- args[-(1:2)]

  if (action == "status") {
    if (length(flags) > 0L) {
      stop(usage, call. = FALSE)
    }
    status <- agent_skill_status()
    if (nrow(status) == 0L) {
      cat("No bixverse agent skill installed. Run `skill install`.\n")
      return(invisible(0L))
    }
    for (i in seq_len(nrow(status))) {
      cat(sprintf(
        "%-6s %-10s %s\n",
        if (status$current[i]) "ok" else "stale",
        if (is.na(status$version[i])) "unknown" else status$version[i],
        status$path[i]
      ))
    }
    cat(sprintf("package: %s\n", utils::packageVersion("bixverse")))
    return(invisible(as.integer(!all(status$current))))
  }

  if (action != "install") {
    stop(usage, call. = FALSE)
  }

  # flags with a value
  value_of <- function(flag) {
    i <- which(flags == flag)
    if (length(i) == 0L) {
      return(NULL)
    }
    if (length(i) > 1L || i == length(flags)) {
      stop(sprintf("`%s` needs exactly one value.", flag), call. = FALSE)
    }
    flags[i + 1L]
  }
  agent <- value_of("--agent")
  dest <- value_of("--dest")
  known <- c("--agent", agent, "--dest", dest, "--project", "--force")
  unknown <- setdiff(flags, known)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown argument(s): %s\n%s",
        paste(unknown, collapse = " "),
        usage
      ),
      call. = FALSE
    )
  }

  install_agent_skill(
    dest = dest,
    agent = if (is.null(agent)) "claude" else agent,
    scope = if ("--project" %in% flags) "project" else "user",
    overwrite = "--force" %in% flags
  )
  invisible(0L)
}
