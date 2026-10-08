"""Desktop scorecards retain exact values and dispatch only supported navigation."""
import json

from starplast import scorecard as S, scorecard_view as V
from starplast.scorecard_browser import ScorecardBrowser


def test_browser_definitions_full_record_export_and_clear():
    card={'scope':{'task':S.T_LABEL,'organism':'synthetic','target':'test'},
        'counts':{'eligible':4,'answered':2,'abstained':2,'correct':1,'wrong':1},
        'metrics':{'accuracy':0.25,'precision_of_calls':0.5,'coverage':0.5}}
    view=V.build_scorecard_view(card,details={'source':{'grade':'synthetic control'},
        'rows':[{'entity':'g1','truth':'a','prediction':None}]},
        links=(V.ScorecardLink('Rows','scorecard:rows/test'),))
    browser=ScorecardBrowser()
    try:
        browser.set_scorecard(view)
        assert 'All eligible: 4' in browser.toPlainText()
        assert 'Calls made: 2' in browser.toPlainText()
        assert browser.navigate('scorecard:metric/accuracy')
        assert 'ALL hidden genes' in browser.toPlainText()
        assert browser.navigate('scorecard:expand')
        assert 'synthetic control' in browser.toPlainText()
        assert json.loads(browser.export_json())['snapshot']['card']==card
        rows=[];browser.rows_requested.connect(rows.append)
        assert browser.navigate('scorecard:rows/test') and rows==['test']
        assert not browser.navigate('javascript:alert(1)')
        assert not browser.navigate('scorecard:detail/unknown')
        assert not browser.navigate('scorecard:rows/unknown')
        assert not browser.navigate('scorecard:outcome/unknown')
        assert not browser.navigate('scorecard:detail/accuracy')
        assert not browser.navigate('scorecard:metric/source')
        assert not browser.navigate('https://unregistered.example/source')
        browser.clear()
        assert browser.view is None and not browser.navigate('scorecard:expand')
    finally:browser.close()
