#!/usr/bin/env python3
"""PreToolUse guard for the physics loop.

Active only in projects that contain notes/STATUS.md. Enforces:
  * tests/acceptance/ is written only by the `planner` agent;
  * each agent writes only to its own folders;
  * notes/STATUS.md is edited only by the main session (coordinator);
  * new files in prompts/ and reports/<kind>/ are named NNN-slug.md with an unused number.
Fails open (allows the call) on any internal error so it can never brick a session.
"""
import json
import os
import re
import shlex
import sys

REPORT_KINDS = {"tests", "benchmarks", "fixes", "audits", "research"}
NAME_RE = re.compile(r"^(\d{3,})-[a-z0-9]+(?:-[a-z0-9]+)*\.md$")
ACC = "tests/acceptance"

# Agents that may write only inside these folders (prefix match, relative to project root).
WRITE_ALLOW = {
    "planner": ["prompts/", "reports/research/", "tests/acceptance/", "docs/", "scratch/"],
    "tester": ["tests/", "reports/tests/", "reports/benchmarks/", "results/"],
    "auditor": ["reports/audits/"],
}
# builder, fixer and any other subagent: anywhere EXCEPT these.
CODE_AGENT_DENY = ["prompts/", "reports/", "notes/STATUS.md", "docs/METHODS.md"]
OWN_REPORTS = {"fixer": "reports/fixes/"}

def deny(reason):
    print(json.dumps({"hookSpecificOutput": {
        "hookEventName": "PreToolUse",
        "permissionDecision": "deny",
        "permissionDecisionReason": "[physics-loop guard] " + reason}}))
    sys.exit(0)

def check_write(rel, role, tool, project, abspath):
    if (rel == ACC or rel.startswith(ACC + "/")) and role != "planner":
        deny("tests/acceptance/ is locked. Only the planner agent may change acceptance tests. "
             "If you think a test is wrong, STOP and report it (see global CLAUDE.md).")
    if role != "main":
        if role in WRITE_ALLOW:
            if not any(rel.startswith(p) for p in WRITE_ALLOW[role]):
                deny(f"The {role} agent may only write under: {', '.join(WRITE_ALLOW[role])} (tried {rel}).")
        else:
            own = OWN_REPORTS.get(role)
            if any(rel.startswith(p) for p in CODE_AGENT_DENY) and not (own and rel.startswith(own)):
                deny(f"The {role} agent may not write {rel}. prompts/, reports/, notes/STATUS.md and "
                     "docs/METHODS.md belong to other roles" + (f"; your reports go in {own}" if own else "") + ".")

    parts = rel.split("/")
    if parts[0] in ("prompts", "reports"):
        if parts[0] == "reports":
            if len(parts) != 3 or parts[1] not in REPORT_KINDS:
                deny("Reports live in reports/<kind>/NNN-slug.md with kind in "
                     f"{sorted(REPORT_KINDS)}. Got: {rel}")
            folder = os.path.join(project, "reports", parts[1])
        else:
            if len(parts) != 2:
                deny(f"Prompts live directly in prompts/ as NNN-slug.md. Got: {rel}")
            folder = os.path.join(project, "prompts")
        frel = os.path.relpath(folder, project)
        m = NAME_RE.match(parts[-1])
        if not m:
            deny(f"Bad name '{parts[-1]}'. Use NNN-short-slug.md (3+ digits, lowercase, hyphens). "
                 f"Get the path with: .claude/bin/next-num {frel} <slug words>")
        if not os.path.exists(abspath) and os.path.isdir(folder):
            used = {int(NAME_RE.match(f).group(1)) for f in os.listdir(folder) if NAME_RE.match(f)}
            if int(m.group(1)) in used:
                deny(f"Number {m.group(1)} is already used in {frel}/. Run .claude/bin/next-num {frel} <slug>")

# ---------- shell commands: block anything that would modify tests/acceptance ----------
WRITERS_DEST = {"cp", "mv", "rsync", "install", "ln"}          # last path argument is written
WRITERS_ANY = {"rm", "rmdir", "truncate", "touch", "chmod", "chown", "unlink", "patch", "tee", "shred"}
GIT_WRITERS = {"checkout", "restore", "rm", "mv", "stash", "apply", "am", "cherry-pick", "revert", "reset", "clean"}

def is_acc(path, project, cwd):
    path = path.strip("'\"")
    if not path or path.startswith("-"):
        return False
    path = os.path.expanduser(path)
    ab = os.path.normpath(path if os.path.isabs(path) else os.path.join(cwd, path))
    rel = os.path.relpath(ab, project).replace(os.sep, "/")
    return rel == ACC or rel.startswith(ACC + "/") or rel.startswith(ACC + "*")

def segments(cmd):
    # split on ; && || | & and newlines, keep redirections inside segments
    return [s for s in re.split(r"\|\||&&|[;|&\n]", cmd) if s.strip()]

def check_bash(cmd, role, project, cwd):
    if role != "main" and "acceptance-lock" in cmd:
        deny("Only the main session (coordinator) runs acceptance-lock.")
    if role == "planner" or "acceptance" not in cmd:
        return
    for redir in re.findall(r"(?:^|[^0-9&<])>{1,2}\s*([^\s;|&<>]+)", cmd):
        if is_acc(redir, project, cwd):
            deny("This command redirects output into tests/acceptance/, which is locked (planner-only).")
    for seg in segments(cmd):
        try:
            toks = shlex.split(seg)
        except ValueError:
            toks = seg.split()
        while toks and (re.match(r"^\w+=", toks[0]) or toks[0] in ("sudo", "env", "command", "nice", "time")):
            toks = toks[1:]
        if not toks:
            continue
        prog = os.path.basename(toks[0])
        args = [t for t in toks[1:] if not t.startswith(">")]
        paths = [a for a in args if not a.startswith("-")]
        hit = False
        if prog in WRITERS_DEST and paths:
            hit = is_acc(paths[-1], project, cwd)
        elif prog in WRITERS_ANY:
            hit = any(is_acc(a, project, cwd) for a in paths)
        elif prog == "dd":
            hit = any(a.startswith("of=") and is_acc(a[3:], project, cwd) for a in args)
        elif prog in ("sed", "perl") and any(re.match(r"^-[a-zA-Z]*i", a) or a.startswith("--in-place") for a in args):
            hit = any(is_acc(a, project, cwd) for a in paths)
        elif prog == "git" and len(toks) > 1 and toks[1] in GIT_WRITERS:
            hit = any(is_acc(a, project, cwd) for a in paths[1:])
        elif prog.startswith("python"):
            hit = bool(re.search(r"open\([^)]*tests/acceptance[^)]*['\"][wax]", seg)
                       or re.search(r"tests/acceptance[^)]*\)\.write_(text|bytes)\(", seg))
        if hit:
            deny("This shell command would modify tests/acceptance/, which is locked. "
                 "Reading or running the tests is fine; changing them is planner-only.")

def main():
    try:
        data = json.load(sys.stdin)
    except Exception:
        sys.exit(0)
    cwd = data.get("cwd") or os.getcwd()
    project = os.path.realpath(os.environ.get("CLAUDE_PROJECT_DIR") or cwd)
    if not os.path.isfile(os.path.join(project, "notes", "STATUS.md")):
        sys.exit(0)  # physics loop not active here
    role = data.get("agent_type") or "main"
    tool = data.get("tool_name", "")
    ti = data.get("tool_input") or {}
    if tool == "Bash":
        check_bash(ti.get("command", ""), role, project, os.path.realpath(cwd))
        sys.exit(0)
    path = ti.get("file_path") or ti.get("notebook_path")
    if not path:
        sys.exit(0)
    abspath = os.path.realpath(path if os.path.isabs(path) else os.path.join(cwd, path))
    rel = os.path.relpath(abspath, project).replace(os.sep, "/")
    if rel.startswith(".."):
        sys.exit(0)  # outside the project: not our business
    check_write(rel, role, tool, project, abspath)
    sys.exit(0)

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception:
        sys.exit(0)
