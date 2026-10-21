#!/usr/bin/env python3
"""SessionStart hook: re-load the physics-loop state into Claude's context.

Runs on startup, resume, /clear and after context compaction. Prints (stdout -> context):
  * notes/STATUS.md (capped), the plan in force, latest reports, acceptance-lock status;
  * refreshes the graphify code graph if one exists.
Silent in folders that are not physics-loop projects.
"""
import os
import re
import subprocess
import sys

NAME_RE = re.compile(r"^(\d{3,})-.*\.md$")

def sh(args, cwd, timeout=15):
    try:
        r = subprocess.run(args, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, (r.stdout + r.stderr).strip()
    except Exception as e:
        return 99, str(e)

def latest(folder, n=1):
    if not os.path.isdir(folder):
        return []
    files = [f for f in os.listdir(folder) if NAME_RE.match(f)]
    return sorted(files, key=lambda f: int(NAME_RE.match(f).group(1)))[-n:]

def main():
    project = os.environ.get("CLAUDE_PROJECT_DIR") or os.getcwd()
    status = os.path.join(project, "notes", "STATUS.md")
    if not os.path.isfile(status):
        if os.path.isdir(os.path.join(project, "prompts")) or os.path.isfile(os.path.join(project, "PLAN.md")):
            print("Note: this folder has a plan (prompts/ or PLAN.md) but the physics loop is not "
                  "initialised (no notes/STATUS.md). If Digonto wants the agent workflow, he can run /kickoff.")
        return
    out = ["=== PHYSICS LOOP ACTIVE (follow the 'physics loop' section of .claude/global-CLAUDE.md) ==="]
    with open(status, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    if len(lines) > 150:
        lines = lines[:150] + [f"... ({len(lines) - 150} more lines in notes/STATUS.md)"]
    out += ["--- notes/STATUS.md ---", *lines, "--- end STATUS ---"]

    p = latest(os.path.join(project, "prompts"))
    out.append(f"Plan in force: prompts/{p[0]}" if p else "Plan in force: none found in prompts/")
    for kind in ("tests", "fixes", "audits", "research", "benchmarks"):
        files = latest(os.path.join(project, "reports", kind), 2)
        if files:
            out.append(f"Latest reports/{kind}: " + ", ".join(files))

    chk = os.path.join(project, ".claude", "bin", "acceptance-check")
    if os.path.isfile(chk):
        code, txt = sh([chk], project)
        out.append("Acceptance lock: " + (txt.splitlines()[0] if txt else f"exit {code}"))

    if os.path.isfile(os.path.join(project, "graphify-out", "graph.json")):
        code, _ = sh(["graphify", "update", "."], project, timeout=12)
        out.append("Code graph: refreshed with `graphify update .`" if code == 0
                   else "Code graph: refresh failed or timed out - run `graphify update .` yourself.")
    else:
        out.append("Code graph: missing - run `graphify update .` before exploring code.")
    print("\n".join(out))

if __name__ == "__main__":
    try:
        main()
    except Exception as e:  # never block a session
        print(f"(physics-loop session hook error: {e})")
    sys.exit(0)
