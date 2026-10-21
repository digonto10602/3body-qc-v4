# graphify
- **graphify** (`~/.claude/skills/graphify/SKILL.md`) - any input to knowledge graph. Trigger: `/graphify`
When the user types `/graphify`, use the installed graphify skill or instructions before doing anything else.

<!-- physics-loop:start -->
# Global rules (Digonto) — loaded in every Claude Code session

## How to report to me (every session, every project)
- Plain English. Define every jargon word the first time you use it (software, tool and physics terms).
- Never invent physics-sounding terms. If a label is your own, say so and explain it.
- Order every report: (1) what was asked, (2) concrete actions taken, (3) what the result means, (4) problems or assumptions.
- All math in LaTeX ($...$ inline, $$...$$ display), never plain-text/unicode math.
- Numbers always come with units, uncertainty/tolerance, and the command or file that produced them.

## Code understanding: graphify first
- If `graphify-out/graph.json` exists, orient with graphify BEFORE grepping or reading many files:
  `graphify query "<question>"`, `graphify explain "<symbol>"`, `graphify path "<A>" "<B>"`,
  `graphify affected "<symbol>"` (what breaks if I change this), `graphify god-nodes` (central code).
- Read raw files only after graphify has oriented you, or to edit/debug specific lines.
- After changing code, run `graphify update .` (local, no API cost) before handing work to another agent.
- Include this rule in every prompt you give a subagent that touches code.

## The physics loop (active when `notes/STATUS.md` exists in the project)
| Role | Who | Model / effort |
|---|---|---|
| coordinator (plans steps, delegates, keeps STATUS, commits) | main session | Opus 5.5, high |
| builder (writes code + unit tests, compiles) | `builder` agent | Sonnet 5.5, high |
| tester (independent verification, benchmarks) | `tester` agent | Sonnet 5.5, high |
| fixer (root-cause debugging) | `fixer` agent | Opus 5.5, high |
| auditor (physics/numerics review at milestones) | `auditor` agent | Opus 5.5, high |
| planner (web research, plans, acceptance tests) | `planner` agent | Fable 5.1, high |

### Coordinator rules (main session)
- Do not write production code yourself beyond trivial one-line edits; delegate to `builder`.
  Your job: pick the next step, brief agents precisely, read their verdicts, update STATUS, commit.
- Brief agents with file paths, test ids, report paths and the plan section — not pasted file contents.
- You are the only one who edits `notes/STATUS.md`, commits, tags, and runs `acceptance-lock`.

### Shared files (the only memory agents share — keep them current)
- `notes/STATUS.md` — current milestone, step, open failures with attempt counters, next action. Update it after every agent result.
- `prompts/NNN-slug.md` — every plan. `prompts/000-*` is the starter plan from chat. Highest number = the plan in force.
- `reports/tests/`, `reports/benchmarks/`, `reports/fixes/`, `reports/audits/`, `reports/research/` — one numbered report per agent run, in its own folder.
- `docs/METHODS.md` — every numerical method in use, its equations, conventions and source (arXiv ID / DOI / URL).
- `results/` — raw numbers (JSON/CSV/HDF5), each with git commit, parameters, seed, machine.
- `tests/acceptance/mN/` — LOCKED definition of "correct" for milestone N. Only the planner writes it; changes are committed + tagged `acceptance-vN` with `~/.claude/bin/acceptance-lock`.

### Numbering rule (enforced by a hook)
- Every new file in `prompts/` or `reports/<kind>/` is named `NNN-short-slug.md` (3+ digits, lowercase, hyphens).
- Get the path with `~/.claude/bin/next-num <folder> <slug words>` (path is relative to the project root). Numbers are per folder, never reused, never renumbered, files never moved between folders.
- Every report starts with: title, date, git commit (`git rev-parse --short HEAD`), the prompt in force, and links to related reports.

### Build loop
1. Read `notes/STATUS.md` and the plan in force. Pick ONE small step with its done-check (the plan lists, per step, which acceptance test ids it must turn green).
2. Call `builder` with: the step, plan section, conventions to respect, files likely involved, the done-check.
3. Call `tester` (mode `step`) with: the step, files changed, test ids to verify.
4. Tester PASS → commit (`git add -A && git commit -m "<step> (prompts/NNN, reports/tests/NNN)"`), update STATUS, next step.
   Tester FAIL → send the failing report to `builder` to fix; record the attempt in STATUS.

### Escalation — objective triggers, counted in STATUS
Call `fixer` when ANY is true (you may escalate earlier if you can see the bug is subtle):
- the same check has failed after 3 builder attempts;
- 2 tester rounds in a row without fewer failing checks;
- the build/compile is still broken after 2 builder attempts;
- a physics check fails and the builder cannot name the cause after 1 attempt (wrong sign, wrong spectrum, broken symmetry, wrong convergence order).
Fixer FIXED → call `tester`; PASS → commit. Call `planner` when ANY is true:
- `fixer` returned NOT FIXED twice on the same issue;
- `fixer` returned NOT FIXED – METHOD (the numerical method or plan is wrong: unstable, not converging, wrong scaling, infeasible cost);
- the step needs a method that `docs/METHODS.md` does not cover.
After a planner run: `~/.claude/bin/acceptance-lock "<plan path>"`, reset the counters in STATUS, continue from step 1 of the new prompt.

### Stop and ask me (do not continue) when
- a planner plan has failed twice (two planner prompts for the same goal both escalated back);
- any agent believes an acceptance test or tolerance itself is wrong;
- a result contradicts established physics and no bug is found;
- you would need to delete data, rewrite git history, push, or spend cluster allocation not mentioned in the prompt.

### Milestone gate (all must pass, in this order)
1. `tester` mode `full`: acceptance tests of milestones 1..N all PASS, unit tests PASS, `~/.claude/bin/acceptance-check` PASS.
2. `auditor` review of the milestone diff → PASS or PASS WITH NOTES (FAIL sends the blocking issues to `fixer`).
3. `tester` mode `benchmark` (only now — never benchmark wrong code); compare with the previous benchmark report.
4. Commit, `git tag milestone-N`, update STATUS, give me a short report (format above).

### Hard rules for correctness
- Never weaken a test, loosen a tolerance, skip, xfail, or special-case a test input to make it pass. Tolerances change only through the planner, with a written justification.
- Never hard-code expected numbers into source code. Never silently catch exceptions in numerical code.
- Fixed random seeds for every stochastic test; record seeds in results.
- Use float64/complex128 for reference computations; any lower precision must be checked against them.
- When CPU and GPU paths both exist, they must agree within the stated tolerance.
- Every number that leaves the project (plot, table, report) must be traceable to a commit + command + parameters.

### Speed rules
- Fast checks after every edit; full suite at milestone gates; benchmarks only after the audit passes.
- Smallest system size that tests the logic. Long jobs go to the background or Slurm, never block the loop.
- Use graphify instead of reading whole files. Do not re-read files already in context.
- If Claude Code offers to continue (resume) the previous builder agent for the next step in the same area, prefer that over a fresh one — it keeps its context.
<!-- physics-loop:end -->
