# 67.01 — Keep the work list running and record session stops

Status: COMPLETE, 2026-10-08; validated setup committed/pushed to nightly.

User request, 2026-10-08: a watchdog that keeps work moving through the list and
logs what makes the session stop. This is separate from the cancelled worker
deployment/retry request.

| Item | Percent done | Time left | Description |
|---|---:|---|---|
| 67.01 | ✅ | 0 h | Enable native goal continuation and independently log activity/session stop reasons |

Bounded implementation: native thread goal controls continuation of instructions
64/65/66. A local once-per-minute systemd timer observes only this thread's log,
with an append-only event log and atomic/locked state. Log actual runtime turn
completions/abort reasons and agent-reported lifecycle stop reasons; flag stale
activity after 15 minutes once, with unknown cause when the runtime supplies none.
Never treat turn completion as whole-list completion. No parallel agent launcher,
database polling, permissions bypass, budget bypass or fabricated heartbeat.

Acceptance: active native goal verified, known-truth stop/stall/partial-write
tests pass, timer installed and active, actual runtime boundary/stop records
verified, durable instructions and logging location documented, commit/push to
nightly. Add this row to the full progress table (49 total); completion requires
all acceptance evidence. Host/bio/release actions keep their existing scope.

Completion evidence: `results/work_watchdog_setup_2026_10_08/`. Native goal
active, timer enabled/active with successful service check, 403 checks passed.
Twelve actual historical boundaries backfilled, four invalid_prompt content
restriction stops recorded. Private live JSONL/state remains ignored by git.
`docs/work-watchdog.md` explains controls, log paths and actual limits.
Instructions 64/65/66 remain active; this completes only the watchdog action.
