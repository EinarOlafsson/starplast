# Continuing the Starplast work list

The native thread goal controls automatic continuation of the unfinished action
list. A separate local timer records runtime boundaries and detects inactivity;
it does not launch another agent. Goals continue at eligible idle boundaries,
subject to user/runtime controls. Official documentation:
https://developers.openai.com/cookbook/examples/codex/using_goals_in_codex

The enabled `starplast-work-watchdog.timer` checks once per minute with a 128 MB
service cap. After 15 minutes without source/progress activity while the last
reported goal state is active, it logs a stale event once. Stale activity cannot
identify a crash, database overload or other unreported cause. Normal turn
completion never implies whole-list completion. Runtime task_complete errors
retain their reported code/message; explicit agent-reported stop/pause/limit
reasons remain separate. Runtime restrictions and explicit interruptions cannot
be overridden by this monitor.

Private files, excluded from git:

- `.starplast-watchdog/events.jsonl`: append-only progress, turn, stop and stale events.
- `.starplast-watchdog/state.json`: current item, source offset, activity and last observed goal status.
- `journalctl --user -u starplast-work-watchdog.service`: service failures, including failures to access state/log files.

At resume, inspect these files and use the goal tool to verify actual lifecycle
state; the monitor's saved goal status is only the latest agent observation.
Record progress and boundaries with explicit reasons:

```sh
python scripts/work_watchdog.py record --kind progress --item 64.11 --goal-status active --reason 'Verified goal active; starting the next frozen numeric adapter'
python scripts/work_watchdog.py record --kind turn_boundary --item 64.11 --goal-status active --reason 'Progress report; unfinished work remains; native continuation expected'
python scripts/work_watchdog.py check
```

Use `--kind stop` for an actual stop with its observed cause/status. Pausing the
native goal requires the user's explicit request; recording a state does not
change the goal. Never write fake periodic heartbeats to hide a stalled agent.
Actual runtime activity is observed independently; no message bodies, reasoning
or tool inputs/outputs are copied. One-time backfill retains only historical
boundary/error metadata. Unknown interruption causes stay unknown.

Installation (already performed on this host):

```sh
cp scripts/systemd/starplast-work-watchdog.service scripts/systemd/starplast-work-watchdog.timer ~/.config/systemd/user/
systemctl --user daemon-reload
systemctl --user enable --now starplast-work-watchdog.timer
systemctl --user start starplast-work-watchdog.service
```

A new thread needs its exact runtime log configured explicitly; an existing
monitor refuses reassignment to another thread. User systemd availability and
native continuation are separate: the timer can record inactivity even when
continuation cannot run. See the acceptance packet at
`results/work_watchdog_setup_2026_10_08/` and completed action 67.01.
