# agent skill ------------------------------------------------------------------

skill_root <- tempfile("skill_test")
dir.create(skill_root)

# tinytest resets the working directory after every expression, so the project
# installs run inside one
in_root <- function(expr) {
  old_wd <- setwd(skill_root)
  on.exit(setwd(old_wd))
  force(expr)
}

pkg_version <- as.character(utils::packageVersion("bixverse"))

## install ---------------------------------------------------------------------

target <- in_root(normalizePath(suppressMessages(install_agent_skill(
  agent = "codex",
  scope = "project"
))))

expect_equal(
  current = normalizePath(target),
  target = normalizePath(file.path(skill_root, ".agents/skills/bixverse")),
  info = "agent skill - codex project install goes to .agents/skills"
)

expect_true(
  current = file.exists(file.path(target, "SKILL.md")),
  info = "agent skill - SKILL.md is kept for codex"
)

expect_equal(
  current = readLines(file.path(target, ".bixverse-version")),
  target = pkg_version,
  info = "agent skill - the install is stamped with the package version"
)

expect_error(
  in_root(install_agent_skill(agent = "codex", scope = "project")),
  pattern = "overwrite",
  info = "agent skill - no silent clobber"
)

## status ----------------------------------------------------------------------

project_row <- function() {
  status <- in_root(agent_skill_status())
  status[in_root(normalizePath(path)) == target]
}

expect_true(
  current = project_row()$current,
  info = "agent skill - a fresh install is current"
)

writeLines("0.0.1", file.path(target, ".bixverse-version"))

expect_false(
  current = project_row()$current,
  info = "agent skill - an older stamp is stale"
)

expect_true(
  current = any(grepl(
    "from version 0.0.1",
    in_root(capture.output(
      bixverse:::.warn_stale_agent_skills(),
      type = "message"
    ))
  )),
  info = "agent skill - attach warns about the stale install"
)

## cli -------------------------------------------------------------------------

expect_equal(
  current = in_root(suppressMessages(skill_cli(c(
    "skill",
    "install",
    "--agent",
    "codex",
    "--project",
    "--force"
  )))),
  target = 0L,
  info = "agent skill cli - install with --force"
)

expect_true(
  current = project_row()$current,
  info = "agent skill cli - the forced install refreshes the stamp"
)

cli_dest <- file.path(skill_root, "elsewhere")
suppressMessages(skill_cli(c("skill", "install", "--dest", cli_dest)))

expect_true(
  current = file.exists(file.path(cli_dest, "bixverse", "SKILL.md")),
  info = "agent skill cli - --dest overrides the agent directory"
)

expect_error(
  skill_cli(c("skill", "install", "--bogus")),
  pattern = "Unknown argument",
  info = "agent skill cli - unknown flags error"
)

expect_error(
  skill_cli(c("skill", "remove")),
  pattern = "Usage",
  info = "agent skill cli - unknown actions error"
)

# clean up ---------------------------------------------------------------------

unlink(skill_root, recursive = TRUE, force = TRUE)
