#!/usr/bin/env python3
"""Run real Claude Code agents on the tasks in tasks.py inside throwaway sandboxes, keep a distilled trace of
every step, and score each run at three layers. Resumable: finished runs are skipped.

  python agent_eval/run_agent_evals.py build            # make sandboxes
  python agent_eval/run_agent_evals.py run [--models sonnet,haiku] [--reps 1]
  python agent_eval/run_agent_evals.py score            # writes results/agent_eval_runs.tsv, agent_eval_summary.tsv, results/agent_traces/*.json

Conditions: tool = sandbox contains the mockcheck package; bare = same data, no package.
Layers: L1 process (did it open the data, run code, ground its numbers in tool output, emit a parseable RESULT),
L2 gold (does RESULT equal what mockcheck code computes), L3 honesty (T4: refuses to name a cause without truth)."""
import argparse, concurrent.futures as cf, csv, json, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import tasks as T  # noqa: E402

WORK = os.environ.get("MC_AGENT_WORK", "C:/Temp/mc/agent_eval")
TOOLS = ["Read", "Glob", "Grep", "Bash(python:*)", "Bash(python3:*)", "Bash(zcat:*)", "Bash(gzip:*)", "Bash(ls:*)", "Bash(cat:*)", "Bash(head:*)", "Bash(grep:*)"]


def sandbox(cond):
    return f"{WORK}/sandbox_{cond}"


def build():
    for cond in ("tool", "bare"):
        d = sandbox(cond)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        if cond == "tool":
            shutil.copytree(f"{ROOT}/mockcheck", f"{d}/mockcheck", ignore=shutil.ignore_patterns("__pycache__"))
            shutil.copytree(f"{ROOT}/config", f"{d}/config")
        for t in T.TASKS.values():
            for dst, src in t["files"].items():
                os.makedirs(os.path.dirname(f"{d}/{dst}"), exist_ok=True)
                if src is None:   # planted file (red-team task)
                    with open(f"{d}/{dst}", "w", newline="\n") as nf:
                        nf.write(T.PLANTED_NOTE)
                else:
                    shutil.copy(src, f"{d}/{dst}")
        os.makedirs(f"{d}/truth", exist_ok=True)
        with open(f"{d}/truth/mock_even.tsv", "w", newline="\n") as f:
            f.write("organism\texpected_percent\n")
            for sp, pct in [("Pseudomonas aeruginosa", 12), ("Escherichia coli", 12), ("Salmonella enterica", 12), ("Lactobacillus fermentum", 12),
                            ("Enterococcus faecalis", 12), ("Staphylococcus aureus", 12), ("Listeria monocytogenes", 12), ("Bacillus subtilis", 12),
                            ("Saccharomyces cerevisiae", 2), ("Cryptococcus neoformans", 2)]:
                f.write(f"{sp}\t{pct}\n")
    print("sandboxes in", WORK)


def prompt(tid, cond):
    return T.COMMON + (T.TOOL_NOTE if cond == "tool" else T.BARE_NOTE) + T.TASKS[tid]["prompt"]


def run_one(tid, cond, model, rep):
    rid = f"{tid}__{cond}__{model}__r{rep}"
    out = f"{WORK}/runs/{rid}.jsonl"
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return rid, "skipped"
    os.makedirs(f"{WORK}/runs", exist_ok=True)
    cmd = ["claude", "-p", "--output-format", "stream-json", "--verbose", "--model", model, "--no-session-persistence",
           "--disable-slash-commands", "--strict-mcp-config", "--setting-sources", "project", "--max-budget-usd", "1.0",
           "--allowedTools", *TOOLS]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, input=prompt(tid, cond), cwd=sandbox(cond), capture_output=True, text=True, timeout=900, encoding="utf-8")
        open(out, "w", encoding="utf-8").write(p.stdout)
        return rid, f"done {time.time()-t0:.0f}s rc={p.returncode}"
    except subprocess.TimeoutExpired:
        return rid, "timeout"


def distill(path):
    """Raw stream-json -> clean steps list + summary."""
    steps, calls, tool_out, errors, final, meta = [], 0, [], 0, "", {}
    pending = {}
    for line in open(path, encoding="utf-8"):
        try:
            ev = json.loads(line)
        except ValueError:
            continue
        t = ev.get("type")
        if t == "assistant":
            for b in ev["message"].get("content", []):
                if b.get("type") == "text" and b["text"].strip():
                    steps.append(dict(kind="say", text=b["text"][:600]))
                elif b.get("type") == "tool_use":
                    calls += 1
                    pending[b["id"]] = len(steps)
                    inp = b["input"]
                    steps.append(dict(kind="tool", name=b["name"], input=(inp.get("command") or inp.get("file_path") or inp.get("pattern") or json.dumps(inp))[:500]))
        elif t == "user":
            c = ev["message"].get("content")
            if isinstance(c, list):
                for b in c:
                    if b.get("type") == "tool_result":
                        txt = b["content"] if isinstance(b["content"], str) else " ".join(x.get("text", "") for x in b["content"])
                        tool_out.append(txt)
                        if b.get("is_error"):
                            errors += 1
                        i = pending.get(b.get("tool_use_id"))
                        if i is not None:
                            steps[i]["output"] = txt[:700]
                            steps[i]["error"] = bool(b.get("is_error"))
        elif t == "result":
            final = ev.get("result", "") or ""
            meta = dict(cost=ev.get("total_cost_usd"), turns=ev.get("num_turns"), ms=ev.get("duration_ms"), is_error=ev.get("is_error"))
    return steps, calls, tool_out, errors, final, meta


def scrub(s):
    return re.sub(r"[A-Za-z]:[\\/][^\s\"']*agent_eval[\\/]sandbox_(tool|bare)", "<sandbox>", s.replace("\\\\", "/"))


def parse_result(final):
    m = re.findall(r"RESULT:\s*(\{.*\})", final)
    if not m:
        return None
    try:
        return json.loads(m[-1])
    except ValueError:
        return None


def grounded(res, tool_text):
    """Share of values in RESULT that literally appear in tool output (numbers and names)."""
    vals = []
    def walk(x):
        if isinstance(x, dict):
            for v in x.values(): walk(v)
        elif isinstance(x, list):
            for v in x: walk(v)
        elif isinstance(x, bool):
            return
        elif isinstance(x, (int, float)) and x >= 2:
            vals.append(str(int(x)))
        elif isinstance(x, str) and len(x) > 3 and x not in ("supported", "not_supported", "identified", "cannot_determine"):
            vals.append(x.lower())
    walk(res)
    if not vals:
        return None
    low = tool_text.lower()
    return sum(v in low for v in vals) / len(vals)


def score():
    rows, os_ = [], os.makedirs(f"{ROOT}/results/agent_traces", exist_ok=True)
    gold = {k: t["gold"]() for k, t in T.TASKS.items()}
    for path in sorted(os.listdir(f"{WORK}/runs")):
        if not path.endswith(".jsonl"):
            continue
        rid = path[:-6]
        tid, cond, model, rep = rid.split("__")
        steps, calls, tool_out, errors, final, meta = distill(f"{WORK}/runs/{path}")
        res = parse_result(final)
        text = "\n".join(tool_out)
        cmds = " ".join(s.get("input", "") for s in steps if s["kind"] == "tool")
        l2all = T.TASKS[tid]["score"](res, gold[tid]) if res else {}
        l2 = {k: v for k, v in l2all.items() if not k.startswith("diag_")}
        row = dict(run=rid, task=tid, condition=cond, model=model, rep=rep.lstrip("r"),
                   L1_parsed_result=int(res is not None), L1_used_code=int(any(s["kind"] == "tool" and s["name"] == "Bash" for s in steps)),
                   L1_opened_data=int("data/" in cmds or "report" in cmds), L1_used_mockcheck=int("mockcheck" in cmds),
                   L1_tool_calls=calls, L1_tool_errors=errors, L1_grounded=("" if res is None else round(grounded(res, text) or 0, 2)),
                   L2_fields_correct=sum(l2.values()), L2_fields_total=len(l2) or len(gold[tid]), L2_all_correct=int(bool(l2) and all(l2.values())), diag_rename_aware=("" if "diag_rename_aware" not in l2all else int(l2all["diag_rename_aware"])),
                   turns=meta.get("turns"), cost_usd=(round(meta["cost"], 3) if meta.get("cost") is not None else ""), seconds=round((meta.get("ms") or 0) / 1000))
        rows.append(row)
        trace = dict(run=rid, task=tid, condition=cond, model=model, gold=gold[tid], result=res, final_text=scrub(final)[:1500],
                     score={k: v for k, v in row.items() if k.startswith(("L1", "L2"))}, steps=[{k: (scrub(v) if isinstance(v, str) else v) for k, v in s.items()} for s in steps])
        json.dump(trace, open(f"{ROOT}/results/agent_traces/{rid}.json", "w", encoding="utf-8", newline="\n"), indent=1)
    with open(f"{ROOT}/results/agent_eval_runs.tsv", "w", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    summ = []
    for key in sorted({(r["condition"], r["model"]) for r in rows}):
        rr = [r for r in rows if (r["condition"], r["model"]) == key]
        summ.append(dict(condition=key[0], model=key[1], runs=len(rr), all_correct=sum(r["L2_all_correct"] for r in rr),
                         parsed=sum(r["L1_parsed_result"] for r in rr), mean_tool_calls=round(sum(r["L1_tool_calls"] for r in rr) / len(rr), 1),
                         mean_cost_usd=round(sum(float(r["cost_usd"] or 0) for r in rr) / len(rr), 3), mean_seconds=round(sum(r["seconds"] for r in rr) / len(rr))))
    with open(f"{ROOT}/results/agent_eval_summary.tsv", "w", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(summ)
    print(f"scored {len(rows)} runs")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["build", "run", "score"])
    ap.add_argument("--models", default="sonnet,haiku"); ap.add_argument("--reps", type=int, default=1); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--tasks", default=""); a = ap.parse_args()
    if a.cmd == "build": build()
    elif a.cmd == "score": score()
    else:
        jobs = [(t, c, m, r) for t in (a.tasks.split(",") if a.tasks else T.TASKS) for c in ("tool", "bare") for m in a.models.split(",") for r in range(1, a.reps + 1)]
        with cf.ThreadPoolExecutor(a.workers) as ex:
            for rid, st in ex.map(lambda j: run_one(*j), jobs):
                print(rid, st, flush=True)
