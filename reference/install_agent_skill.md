# Install the bixverse agent skill

bixverse ships a skill that teaches coding agents how to *use* the
package: the workflows, the `params_*()` convention, the ordering
constraints and the traps. This copies it out of the installed package
to where your agent looks for it.

The skill is versioned with the package. The install is stamped with the
package version, and
[`library(bixverse)`](https://gregorlueg.github.io/bixverse) warns when
an installed skill is older than the package; re-run this with
`overwrite = TRUE` after upgrading. From a shell, see
[`skill_cli()`](https://gregorlueg.github.io/bixverse/reference/skill_cli.md).

## Usage

``` r
install_agent_skill(
  dest = NULL,
  agent = c("claude", "codex"),
  scope = c("user", "project"),
  overwrite = FALSE
)
```

## Arguments

- dest:

  String or `NULL`. Skills directory to install into; the skill lands in
  a `bixverse` subdirectory. Overrides `agent` and `scope`.

- agent:

  String. One of `c("claude", "codex")`.

- scope:

  String. One of `c("user", "project")`: under the home directory, or
  under the working directory.

- overwrite:

  Boolean. Replace an existing installation. Defaults to `FALSE`, in
  which case an existing `bixverse` directory causes an error rather
  than a silent clobber.

## Value

The path the skill was written to, invisibly.

## Details

`SKILL.md` follows the Agent Skills open standard, so the same files
serve every agent; only the directory differs. `"claude"` installs to
`.claude/skills`, `"codex"` to `.agents/skills`, which Codex and other
agents scan. A project install ties the skill to the project, so two
projects on different bixverse versions each see their own API. Run it
from the project root.

## Examples

``` r
if (FALSE) { # \dontrun{
# Claude Code, into ~/.claude/skills
install_agent_skill()

# Codex and other open standard agents, into ./.agents/skills
install_agent_skill(agent = "codex", scope = "project")

# refresh after a package upgrade
install_agent_skill(overwrite = TRUE)
} # }
```
