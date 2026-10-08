"""Monitor one Starplast thread without launching workers or changing its runtime."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STATE = ROOT / '.starplast-watchdog'
STATUSES = {'active', 'paused', 'complete', 'blocked', 'budget_limited', 'usage_limited'}
ENDS = {'task_complete': 'turn_completed', 'turn_aborted': 'runtime_turn_aborted',
        'task_failed': 'runtime_task_failed', 'session_end': 'runtime_session_end'}


def _time():
    return datetime.now(timezone.utc).timestamp()


def _stamp(value):
    return datetime.fromtimestamp(value, timezone.utc).isoformat()


@contextmanager
def _locked(folder):
    folder.mkdir(parents=True, exist_ok=True)
    with (folder/'lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def _save(folder, state):
    # Separate writers (timer and agent) serialize through the same file lock.
    descriptor, name = tempfile.mkstemp(dir=folder, prefix='state-')
    try:
        with os.fdopen(descriptor, 'w') as file:
            json.dump(state, file, indent=2, allow_nan=False)
            file.write('\n'); file.flush(); os.fsync(file.fileno())
        os.replace(name, folder/'state.json')
    finally:
        if Path(name).exists():
            Path(name).unlink()


def _event(folder, state, kind, now, **detail):
    record = {'timestamp': _stamp(now), 'thread_id': state['thread_id'],
              'event': kind, 'goal_status_last_reported': state['goal_status'], **detail}
    with (folder/'events.jsonl').open('a') as file:
        file.write(json.dumps(record, allow_nan=False)+'\n')
        file.flush(); os.fsync(file.fileno())


def _boundary(folder, state, entry, now, *, historical=False):
    payload = entry.get('payload', {})
    if entry.get('type') != 'event_msg' or not isinstance(payload, dict) or payload.get('type') not in ENDS:
        return False
    kind = payload['type']
    reason = payload.get('reason')
    error = payload.get('error')
    error = error if isinstance(error, dict) else {}
    code = error.get('codex_error_info')
    message = error.get('message')
    _event(folder, state, 'runtime_boundary', now, source_timestamp=entry.get('timestamp'),
        source_event=kind, reason='runtime_reported_error' if error else ENDS[kind],
        reported_reason=reason if isinstance(reason, str) else None,
        reported_error_code=code if isinstance(code, (str, dict)) else None,
        reported_error_message=message[:2048] if isinstance(message, str) else None,
        historical=historical, interpretation='turn completion is not list/goal completion')
    return True


def backfill(folder, *, now=None):
    """Copy only known historical boundary/error metadata once; omit transcript bodies."""
    now = _time() if now is None else now
    folder = Path(folder)
    with _locked(folder):
        state = json.loads((folder/'state.json').read_text())
        if state.get('history_backfilled'):
            return {'already_backfilled': True}
        count = 0
        with Path(state['source']).open('rb') as source:
            # Capture only the history preceding configuration, avoiding overlap
            # with the live monitor and retaining bounded memory per line.
            remaining = state['offset_at_configuration']
            discard = False
            while remaining:
                line = source.readline(min(remaining, 4*1024*1024))
                if not line:
                    break
                remaining -= len(line)
                complete = line.endswith(b'\n')
                if discard:
                    discard = not complete
                    continue
                if not complete:
                    discard = True
                    continue
                try:
                    entry = json.loads(line)
                except (ValueError, TypeError):
                    continue
                count += int(_boundary(folder, state, entry, now, historical=True))
        state['history_backfilled'] = True
        _event(folder, state, 'history_backfilled', now, boundary_records=count,
            reason='Known historical runtime boundary/error metadata copied once; transcript bodies omitted')
        _save(folder, state)
        return {'historical_boundary_records': count}


def configure(folder, thread_id, source, *, now=None, stale_seconds=900):
    """Bind one existing runtime log; never replace an active monitor silently."""
    now = _time() if now is None else now
    folder, source = Path(folder), Path(source).resolve()
    if stale_seconds < 60 or not thread_id or not source.is_file() or thread_id not in source.name:
        raise ValueError('Provide the matching thread log and a threshold of at least 60 seconds')
    with _locked(folder):
        if (folder/'state.json').exists():
            old = json.loads((folder/'state.json').read_text())
            if old['thread_id'] != thread_id or old['source'] != str(source):
                raise ValueError('Existing watchdog belongs to a different session')
            return old
        stat = source.stat()
        state = {'schema_version': 1, 'thread_id': thread_id, 'source': str(source),
            'source_inode': stat.st_ino, 'offset': stat.st_size, 'offset_at_configuration': stat.st_size, 'goal_status': 'active',
            'last_activity': now, 'stale_seconds': stale_seconds, 'alerted': False,
            'current_item': '67.01', 'latest_reason': 'watchdog_configured'}
        _event(folder, state, 'configured', now, reason='native_goal_active; runtime log monitor starts at current EOF')
        _save(folder, state)
        return state


def record(folder, kind, reason, *, item='', goal_status=None, now=None):
    """Record observed progress/boundaries/stops with an explicit reason supplied by the agent."""
    if kind not in {'progress', 'turn_boundary', 'stop', 'goal_status'} or not reason.strip():
        raise ValueError('Supply an explicit supported event and reason')
    if goal_status is not None and goal_status not in STATUSES:
        raise ValueError('Unknown goal status')
    now = _time() if now is None else now
    folder = Path(folder)
    with _locked(folder):
        state = json.loads((folder/'state.json').read_text())
        if goal_status is not None:
            state['goal_status'] = goal_status
        state['last_activity'] = now
        state['latest_reason'] = reason
        if item:
            state['current_item'] = item
        state['alerted'] = False
        _event(folder, state, kind, now, reason=reason, item=state['current_item'])
        _save(folder, state)
        return state


def check(folder, *, now=None):
    """Consume bounded runtime events and log stale activity once, without inventing a cause."""
    now = _time() if now is None else now
    folder = Path(folder)
    with _locked(folder):
        state = json.loads((folder/'state.json').read_text())
        source = Path(state['source'])
        if not source.exists():
            if not state.get('source_missing'):
                _event(folder, state, 'source_missing', now, reason='runtime event log unavailable; cause unknown')
            state['source_missing'] = True
        else:
            state['source_missing'] = False
            stat = source.stat()
            if stat.st_ino != state['source_inode'] or stat.st_size < state['offset']:
                _event(folder, state, 'source_reset', now, reason='runtime log replaced or truncated; not evidence of task completion')
                state.update(offset=0, source_inode=stat.st_ino)
            with source.open('rb') as file:
                file.seek(state['offset'])
                chunk = file.read(4*1024*1024)
            if stat.st_mtime <= now:
                state['last_activity'] = max(state['last_activity'], stat.st_mtime)
            if state.get('discarding_oversized_line'):
                boundary = chunk.find(b'\n')
                if boundary < 0:
                    state['offset'] += len(chunk)
                    chunk = b''
                else:
                    state['offset'] += boundary+1
                    chunk = chunk[boundary+1:]
                    state['discarding_oversized_line'] = False
            # Never consume a partially written final JSON line.
            length = chunk.rfind(b'\n')+1
            if not length and len(chunk) == 4*1024*1024:
                _event(folder, state, 'oversized_runtime_record', now,
                    reason='Record exceeds bounded read size; skip content, retain source activity and subsequent boundaries')
                state['offset'] += len(chunk)
                state['discarding_oversized_line'] = True
            for line in chunk[:length].splitlines():
                try:
                    entry = json.loads(line)
                    timestamp = datetime.fromisoformat(entry['timestamp'].replace('Z', '+00:00')).timestamp()
                except (ValueError, KeyError, TypeError):
                    _event(folder, state, 'unparsed_runtime_record', now, reason='record schema/timestamp invalid; content not copied')
                    continue
                state['last_activity'] = max(state['last_activity'], min(now, timestamp))
                _boundary(folder, state, entry, now)
            state['offset'] += length
        age = max(0., now-state['last_activity'])
        stale = age >= state['stale_seconds'] and state['goal_status'] == 'active'
        if stale and not state['alerted']:
            _event(folder, state, 'activity_stale', now, seconds_without_activity=age,
                item=state['current_item'], reason='No fresh runtime/progress activity; stop cause unknown',
                continuation='Native active goal controls resumption; monitor never starts a competing agent')
        if not stale and state['alerted']:
            _event(folder, state, 'activity_resumed', now, reason='Fresh activity observed or goal no longer active')
        state['alerted'] = stale
        _save(folder, state)
        return {'goal_status': state['goal_status'], 'current_item': state['current_item'],
                'stale': stale, 'seconds_without_activity': age, 'source_missing': state['source_missing']}


def main():
    """Configure, check or annotate the local session monitor."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state-dir', type=Path, default=DEFAULT_STATE)
    commands = parser.add_subparsers(dest='command', required=True)
    setup = commands.add_parser('configure')
    setup.add_argument('--thread', required=True)
    setup.add_argument('--source', type=Path, required=True)
    setup.add_argument('--stale-seconds', type=int, default=900)
    commands.add_parser('check')
    commands.add_parser('backfill')
    note = commands.add_parser('record')
    note.add_argument('--kind', choices=('progress','turn_boundary','stop','goal_status'), required=True)
    note.add_argument('--reason', required=True)
    note.add_argument('--item', default='')
    note.add_argument('--goal-status', choices=sorted(STATUSES))
    args = parser.parse_args()
    if args.command == 'configure':
        result = configure(args.state_dir,args.thread,args.source,stale_seconds=args.stale_seconds)
    elif args.command == 'check':
        result = check(args.state_dir)
    elif args.command == 'backfill':
        result = backfill(args.state_dir)
    else:
        result = record(args.state_dir,args.kind,args.reason,item=args.item,goal_status=args.goal_status)
    print(json.dumps(result,allow_nan=False))


if __name__ == '__main__':
    main()
