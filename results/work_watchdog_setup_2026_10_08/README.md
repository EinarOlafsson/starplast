# Work continuation watchdog acceptance, 2026-10-08

The user explicitly requested continued work and a watchdog with stop-cause
logging. Native goal continuation is enabled for the unfinished Starplast list.
The native goal is observed active; it is not marked complete by this setup.
Official behavior: https://developers.openai.com/cookbook/examples/codex/using_goals_in_codex

Independent systemd user timer `starplast-work-watchdog.timer` is enabled/active,
checks every minute; its oneshot service succeeds under MemoryMax=128M and
TimeoutStartSec=30. Installed unit copies match `scripts/systemd/`. The monitor
is observational: it never starts a second agent, changes budgets/permissions,
writes the source runtime log or sends prompts to evade a restriction.

Private runtime state/log: `.starplast-watchdog/state.json` and `events.jsonl`,
ignored by git. Locked append-only event writes and atomic state replacement
serialize timer/agent updates. Bounded reads preserve incomplete lines and skip
oversized bodies without losing subsequent boundaries. Real source activity
clears stale observations; 15-minute inactivity emits one event with unknown
cause rather than claiming crash, overload or completion. Explicit reported
pause/limit/completion states are retained and not counted as stalls.

The initial backfill records 12 actual historical runtime boundaries, including
four `invalid_prompt` errors reporting content restrictions. Only boundary/error
metadata is copied, never private agent/user message bodies, reasoning or tool
inputs/outputs. A `task_complete` carrying an error is classified as a reported
runtime failure; a normal turn ending is not goal/list completion. Runtime
failure logging cannot override the failure or ensure resumption after an
explicit interruption, unavailable connection or account/runtime limit.

403 checks passed (watchdog truth controls, public docstrings and organism
invariants), plus systemd unit verification and real service/log inspection.
`verification.json` records the acceptance observations. Code/unit/tests hashes
are frozen in `packet_sha256.json`; live private logs remain outside this packet.
