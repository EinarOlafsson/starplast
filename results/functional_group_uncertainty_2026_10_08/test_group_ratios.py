"""Group resampling retains paired/conditional denominators and exact replay."""
import numpy as np
import pytest

from starplast import group_ratios as G


def fixture():
    return (['a', 'a', 'a', 'b', 'c', 'd', 'e'],
            [[1, 1], [0, 0], [1, 1], [0, 0], [0, 0], [0, 0], [1, 1]],
            [[1, 1], [1, 0], [1, 1], [1, 0], [1, 0], [1, 0], [1, 1]], ['all', 'conditional'])


def test_whole_groups_share_draws_and_replay_exactly():
    state, arrays = G.estimate(*fixture(), seed=23, draws=50)
    groups, a, b, names = fixture()
    for i, selected in enumerate(arrays['selections']):
        positions = [j for g in selected for j, name in enumerate(groups) if name == state['group_order'][g]]
        for column in range(2):
            numerator = sum(a[j][column] for j in positions)
            denominator = sum(b[j][column] for j in positions)
            assert arrays['numerators'][i, column] == numerator
            assert arrays['denominators'][i, column] == denominator
    replay, values = G.estimate(*fixture(), seed=23, draws=50, selections=arrays['selections'])
    assert state == replay
    for name in arrays: np.testing.assert_array_equal(arrays[name], values[name])
    for i, name in enumerate(names):
        valid = arrays['ratios'][:, i][arrays['denominators'][:, i] != 0]
        assert state['rates'][name]['interval'] == np.quantile(valid, [.025, .975], method='linear').tolist()
    assert state['rates']['all']['value'] == 3 / 7
    assert state['rates']['conditional']['value'] == 1
    assert state['rates']['conditional']['undefined_draws'] > 0


def test_zero_denominator_small_groups_and_signed_paired_difference():
    groups = list('abcde')
    state, arrays = G.estimate(groups, [[0, -1]] * 5, [[0, 1]] * 5, ['none', 'difference'], seed=1)
    assert state['rates']['none']['interval'] is None
    assert state['rates']['none']['value'] is None
    assert state['rates']['none']['undefined_draws'] == 500
    assert state['rates']['difference']['interval'] == [-1, -1]
    state, arrays = G.estimate(['one'], [[1]], [[1]], ['rate'], seed=1)
    assert state['status'] == 'unavailable' and state['draws_evaluated'] == 0
    assert state['rates']['rate']['interval'] is None
    assert arrays['selections'].shape == (0, 1)


@pytest.mark.parametrize('mutation', ['group', 'shape', 'float', 'negative_denominator', 'unbounded', 'duplicate_name', 'seed', 'draws', 'replay'])
def test_unresolved_or_invalid_inputs_are_refused(mutation):
    groups, a, b, names = fixture(); kwargs = dict(seed=23, draws=50)
    if mutation == 'group': groups[0] = None
    elif mutation == 'shape': a.pop()
    elif mutation == 'float': a[0][0] = .5
    elif mutation == 'negative_denominator': b[0][0] = -1
    elif mutation == 'unbounded': a[0][0] = 2
    elif mutation == 'duplicate_name': names[1] = names[0]
    elif mutation == 'seed': kwargs['seed'] = True
    elif mutation == 'draws': kwargs['draws'] = 1
    else: kwargs['selections'] = np.full((50, 5), 5)
    with pytest.raises(ValueError): G.estimate(groups, a, b, names, **kwargs)
