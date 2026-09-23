"""Paired uncertainty comparisons must preserve evaluation membership and truth."""
import pandas as pd
import pytest

from scripts.summarize_inference import paired_interval


def evaluation():
    """Two independent groups with two outcomes each."""
    return pd.DataFrame({
        'gene_id': ['a', 'b', 'c', 'd'], 'group': ['first', 'first', 'second', 'second'],
        'fold': [0, 0, 1, 1], 'truth': ['x', 'y', 'x', 'y'],
        'prediction': ['x', 'x', 'y', 'y'], 'role': ['held_out_evaluation'] * 4,
    })


def test_identical_predictions_have_zero_paired_uncertainty():
    frame = evaluation()
    result = paired_interval(frame, frame.sample(frac=1, random_state=4), True, repetitions=20)
    assert result['difference'] == result['low_95'] == result['high_95'] == 0


@pytest.mark.parametrize('column,value', [('group', 'changed'), ('fold', 3), ('truth', 'y')])
def test_unmatched_evaluations_cannot_be_presented_as_paired(column, value):
    frame = evaluation()
    other = frame.copy()
    other.loc[0, column] = value
    with pytest.raises(ValueError, match='paired comparison'):
        paired_interval(frame, other, True, repetitions=20)


def test_regression_difference_has_candidate_minus_baseline_direction():
    baseline = evaluation().assign(truth=[0., 0., 0., 0.], prediction=[2., 2., 2., 2.])
    candidate = baseline.assign(prediction=1.)
    result = paired_interval(baseline, candidate, False, repetitions=20)
    assert result['difference'] == result['low_95'] == result['high_95'] == -1
