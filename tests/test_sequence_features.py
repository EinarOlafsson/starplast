"""Long-protein encoding preserves every residue exactly once in pooled weight."""
import numpy as np
import pytest
from starplast.sequence_features import windows, cache_key


@pytest.mark.parametrize('length',[1,1000,1001,2100,12000])
def test_window_weights_cover_each_residue_once(length):
    coverage=np.zeros(length)
    for a,b,w in windows(length):
        assert 0 <= a < b <= length and b-a <= 1000
        coverage[a:b]+=w
    np.testing.assert_allclose(coverage,1)


def test_nonoverlapping_windows_preserve_terminal_fragment():
    parts=windows(2100,1000,0)
    assert [(a,b) for a,b,_ in parts]==[(0,1000),(1000,2000),(2000,2100)]


@pytest.mark.parametrize('length,size,overlap',[(0,1000,128),(10,1000,1000),(10,1023,128)])
def test_invalid_window_contract_is_refused(length,size,overlap):
    with pytest.raises(ValueError): windows(length,size,overlap)


def test_cache_identity_includes_sequence_revision_and_pooling():
    keys={cache_key('ACDE'),cache_key('ACDF'),cache_key('ACDE',revision='different'),
          cache_key('ACDE',overlap=0),cache_key('ACDE',size=512)}
    assert len(keys)==5
