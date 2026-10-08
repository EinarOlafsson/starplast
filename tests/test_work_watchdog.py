"""Known runtime boundaries, stale detection and explicit stop-cause preservation."""
from datetime import datetime, timezone
import json

import pytest

from scripts import work_watchdog as W


def setup(tmp_path):
    source = tmp_path/'rollout-thread-test.jsonl'
    source.write_text('')
    folder = tmp_path/'state'
    W.configure(folder,'thread-test',source,now=1000,stale_seconds=60)
    return folder,source


def append(source, time, kind, **payload):
    with source.open('a') as file:
        file.write(json.dumps({'timestamp':datetime.fromtimestamp(time,timezone.utc).isoformat(),
            'type':'event_msg','payload':{'type':kind,**payload}})+'\n')


def events(folder):
    return [json.loads(line) for line in (folder/'events.jsonl').read_text().splitlines()]


def test_real_turn_completion_is_logged_without_marking_goal_complete(tmp_path):
    folder,source=setup(tmp_path)
    append(source,1010,'task_complete',last_assistant_message='Private content must not be copied')
    assert not W.check(folder,now=1020)['stale']
    boundary=events(folder)[-1]
    assert boundary['source_event']=='task_complete' and boundary['reason']=='turn_completed'
    assert boundary['goal_status_last_reported']=='active'
    assert 'Private content' not in (folder/'events.jsonl').read_text()


def test_unknown_inactivity_is_deduplicated_and_actual_activity_clears_it(tmp_path):
    folder,source=setup(tmp_path)
    assert W.check(folder,now=1060)['stale']
    assert W.check(folder,now=1070)['stale']
    assert sum(e['event']=='activity_stale' for e in events(folder))==1
    assert events(folder)[-1]['reason'].endswith('cause unknown')
    append(source,1080,'token_count')
    assert not W.check(folder,now=1081)['stale']
    assert events(folder)[-1]['event']=='activity_resumed'


def test_reported_abort_reason_survives_and_is_not_inferred_from_missing_activity(tmp_path):
    folder,source=setup(tmp_path)
    append(source,1010,'turn_aborted',reason='user_interrupt')
    W.check(folder,now=1011)
    assert events(folder)[-1]['reported_reason']=='user_interrupt'
    assert events(folder)[-1]['reason']=='runtime_turn_aborted'


@pytest.mark.parametrize('status',['paused','complete','blocked','budget_limited','usage_limited'])
def test_explicit_lifecycle_stop_is_preserved_and_not_treated_as_a_stall(tmp_path,status):
    folder,_=setup(tmp_path)
    W.record(folder,'stop','Observed lifecycle '+status,goal_status=status,item='64.11',now=1001)
    assert not W.check(folder,now=9000)['stale']
    assert events(folder)[-1]['reason']=='Observed lifecycle '+status


def test_partial_json_is_not_consumed_and_later_boundary_is_read_once(tmp_path):
    folder,source=setup(tmp_path)
    line=json.dumps({'timestamp':datetime.fromtimestamp(1010,timezone.utc).isoformat(),
        'type':'event_msg','payload':{'type':'task_complete'}})
    source.write_text(line[:30])
    W.check(folder,now=1011)
    assert json.loads((folder/'state.json').read_text())['offset']==0
    with source.open('a') as file:file.write(line[30:]+'\n')
    W.check(folder,now=1012);W.check(folder,now=1013)
    assert sum(e['event']=='runtime_boundary' for e in events(folder))==1


def test_missing_source_and_truncation_are_recorded_without_claiming_a_cause(tmp_path):
    folder,source=setup(tmp_path)
    append(source,1010,'token_count');W.check(folder,now=1011)
    source.write_text('')
    W.check(folder,now=1012)
    assert events(folder)[-1]['event']=='source_reset'
    source.unlink()
    assert W.check(folder,now=1013)['source_missing']
    W.check(folder,now=1014)
    assert sum(e['event']=='source_missing' for e in events(folder))==1


def test_a_different_thread_cannot_reconfigure_the_watchdog(tmp_path):
    folder,_=setup(tmp_path)
    other=tmp_path/'rollout-thread-other.jsonl';other.write_text('')
    with pytest.raises(ValueError,match='different session'):
        W.configure(folder,'thread-other',other,now=1020)


def test_oversized_runtime_content_cannot_hide_later_stop_boundaries(tmp_path):
    folder,source=setup(tmp_path)
    source.write_text('x'*(4*1024*1024+100)+'\n')
    append(source,1010,'turn_aborted',reason='reported_test_interrupt')
    W.check(folder,now=1011)
    assert events(folder)[-1]['event']=='oversized_runtime_record'
    W.check(folder,now=1012)
    assert events(folder)[-1]['reported_reason']=='reported_test_interrupt'


def test_task_complete_with_a_reported_error_is_a_failure_not_normal_completion(tmp_path):
    folder,source=setup(tmp_path)
    append(source,1010,'task_complete',error={'message':'Runtime restriction', 'codex_error_info':'invalid_prompt'})
    W.check(folder,now=1011)
    boundary=events(folder)[-1]
    assert boundary['reason']=='runtime_reported_error'
    assert boundary['reported_error_code']=='invalid_prompt'
    assert boundary['reported_error_message']=='Runtime restriction'
    assert boundary['goal_status_last_reported']=='active'


def test_backfill_records_only_history_and_never_changes_live_activity(tmp_path):
    source=tmp_path/'rollout-thread-test.jsonl';source.write_text('')
    append(source,500,'task_complete',last_agent_message='Private history omitted')
    folder=tmp_path/'state'
    W.configure(folder,'thread-test',source,now=1000,stale_seconds=60)
    append(source,1010,'turn_aborted',reason='live_interrupt')
    assert W.backfill(folder,now=1011)['historical_boundary_records']==1
    assert W.backfill(folder,now=1012)['already_backfilled']
    assert json.loads((folder/'state.json').read_text())['last_activity']==1000
    assert 'Private history' not in (folder/'events.jsonl').read_text()
    W.check(folder,now=1013)
    boundaries=[e for e in events(folder) if e['event']=='runtime_boundary']
    assert len(boundaries)==2 and boundaries[0]['historical'] and not boundaries[1]['historical']
