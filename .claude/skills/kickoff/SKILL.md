---
name: kickoff
description: Start the physics loop in the current project - scaffold folders, build the graphify code graph, have the planner turn the starter plan into a numbered plan plus locked acceptance tests, then begin building. Run once per project.
argument-hint: "[path to starter plan] [--review]"
disable-model-invocation: true
---
# /kickoff — start the physics loop here

Arguments: `$ARGUMENTS` (optional path to the starter plan from chat; `--review` = pause for
Digonto's approval before building). You are the coordinator (main session). Follow these steps
in order and report in Digonto's format (what was asked → actions → what it means → problems/assumptions).

## 1. Scaffold (idempotent — never overwrite existing files)
```bash
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || git init
mkdir -p prompts reports/tests reports/benchmarks reports/fixes reports/audits reports/research \
         notes docs tests/acceptance tests/unit results scratch
T=~/.claude/templates/physics-project   # not in this repo; every file it would copy already exists here
[ -f notes/STATUS.md ]  || cp "$T/notes/STATUS.md" notes/STATUS.md
[ -f docs/METHODS.md ]  || cp "$T/docs/METHODS.md" docs/METHODS.md
[ -f CLAUDE.md ]        || cp "$T/CLAUDE.md" CLAUDE.md
for p in graphify-out/ scratch/ __pycache__/ '*.pyc' .pytest_cache/; do
  grep -qxF "$p" .gitignore 2>/dev/null || echo "$p" >> .gitignore
done
```

## 2. Starter plan → prompts/000
- Find the starter plan: the argument path, else an existing `prompts/000-*.md`, else `PLAN.md`,
  else ask Digonto.
- If no `prompts/000-*.md` exists, copy (do not move) the starter plan to
  `prompts/000-starter-plan.md`. Other unnumbered plan files from the starter package get the next
  numbers from `.claude/bin/next-num prompts <slug>` (copy, keep the originals).
- Fill the project `CLAUDE.md` placeholders (goal, environment, commands, layout) from the
  starter package. Keep it short. The fast-check command must exclude `tests/acceptance`.

## 3. Code graph
```bash
graphify update .        # local AST graph, no API cost (skip with a warning if graphify is missing)
graphify hook install    # rebuild the graph automatically on every commit / checkout
```

## 4. First commit
`git add -A && git commit -m "kickoff: physics-loop scaffolding + starter plan"`

## 5. Planner (kickoff mode)
Call the `planner` agent with: "KICKOFF. Starter plan: <the actual prompts/000-… file>. Check the
plan's physics and feasibility, research the best numerical methods, and write the first numbered
plan (milestones, small steps with done-checks and acceptance test ids, Interface section), the
acceptance tests in tests/acceptance/mN/, and docs/METHODS.md. Use graphify for any existing code."

## 6. Stubs and lock
- Call the `builder` agent: "Create the interface stubs listed in the Interface section of
  <plan path> (functions/classes that raise NotImplementedError). No other code."
- Run the acceptance tests yourself: every test must FAIL with `NotImplementedError` — not error
  on import or syntax. If not, send the exact output to the planner once; if still broken, stop and ask.
- `.claude/bin/acceptance-lock "kickoff: <plan path>"`

## 7. STATUS and go
- Fill `notes/STATUS.md`: plan in force, acceptance tag, Milestone 1 / Step 1, next action.
- `git add -A && git commit -m "kickoff: plan, stubs, locked acceptance tests"`
- Report to Digonto: the method chosen (with sources), milestones, what the acceptance tests
  check, anything that needs his decision.
- If `--review` was given, STOP here and wait. Otherwise start the build loop at Step 1.
