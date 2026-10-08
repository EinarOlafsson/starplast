"""Native ridge parity, preprocessing leakage and unsupported cohorts."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from starplast import organisms as O, strategies as S, strategy_catalog as N, value_records as V
from starplast.splits import Assignment, ExclusionManifest, SplitManifest

IDS = tuple(f'TGME49_{i:06d}' for i in range(100001, 100036))
ROLES = ['train']*24 + ['tune']*3 + ['calibration']*3 + ['test']*5
SPLIT = SplitManifest(O.TOXOPLASMA, 'ridge-fixture', 'entity', tuple(Assignment(e, e, r) for e, r in zip(IDS, ROLES)), 17)
EXCLUSIONS = ExclusionManifest(SPLIT.benchmark_id, ('target',), ('target', 'sibling'), (), SPLIT.entities('train'), SPLIT.identity)


def inputs():
    features = pd.DataFrame({'x': np.arange(35.)%6, 'y': np.linspace(0, 2, 35), 'empty': np.nan}, index=IDS)
    features.loc[IDS[2], 'x'] = np.nan
    features.loc[IDS[-1], 'y'] = 4.
    values = pd.Series(np.sin(np.arange(24.)), index=SPLIT.entities('train'))
    return features, values


def test_training_matrix_and_predictions_match_native_with_frozen_query_ranks():
    features, values = inputs()
    batch = V.trait_ridge(features, values, split=SPLIT, exclusions=EXCLUSIONS)
    state = batch.model_state
    native_context = S.Context(features.loc[list(values.index), ['x','y']].reset_index(drop=True).assign(gene_id=list(values.index)), graph={}, organism=O.TOXOPLASMA)
    np.testing.assert_array_equal(state['training_vectors'], native_context.matrix(['x','y']))
    # Independent application of the saved rank distributions.
    matrix = np.column_stack([np.array([((np.searchsorted(state['training_distributions'][c],v,'left')+
        np.searchsorted(state['training_distributions'][c],v,'right')+1)/2 if v in state['training_distributions'][c]
        else np.searchsorted(state['training_distributions'][c],v,'right'))/len(state['training_distributions'][c])-.5
        if np.isfinite(v) else 0. for v in features[c]]) for c in state['kept_columns']])
    visible = pd.Series([*values, *([np.nan]*11)])
    native = N._fit_predict(matrix, visible, 'ridge')
    np.testing.assert_array_equal(batch.rows.prediction, native[-5:])
    assert batch.rows.entity.tolist()==list(SPLIT.entities('test'))
    assert not batch.rows.abstained.any() and batch.rows.calibrated_confidence.isna().all()
    assert state['withheld_all_missing_training_columns']==['empty']


def test_held_out_features_cannot_change_fitted_ranks_or_coefficients():
    features, values = inputs()
    original = V.trait_ridge(features, values, split=SPLIT, exclusions=EXCLUSIONS)
    changed = features.copy(); changed.loc[list(SPLIT.entities('test')), 'x'] = 1e8
    altered = V.trait_ridge(changed, values, split=SPLIT, exclusions=EXCLUSIONS)
    assert altered.model_state==original.model_state
    assert not np.array_equal(altered.rows.prediction, original.rows.prediction)


def test_target_family_or_nontraining_values_are_refused():
    features, values = inputs()
    with pytest.raises(ValueError, match='contamination'):
        V.trait_ridge(features.assign(target=0.), values, split=SPLIT, exclusions=EXCLUSIONS)
    with pytest.raises(ValueError, match='training values only'):
        V.trait_ridge(features, values.iloc[::-1], split=SPLIT, exclusions=EXCLUSIONS)
    with pytest.raises(ValueError, match='ordered training cohort'):
        V.trait_ridge(features, values, split=SPLIT, exclusions=replace(EXCLUSIONS,fit_entities=SPLIT.entities('test')))
    with pytest.raises(ValueError, match='population/order'):
        V.trait_ridge(features.iloc[::-1], values, split=SPLIT, exclusions=EXCLUSIONS)


def test_invalid_numbers_and_duplicate_inputs_are_refused():
    features, values = inputs()
    for value in (np.nan, np.inf):
        bad = values.copy(); bad.iloc[0]=value
        with pytest.raises(ValueError, match='observed finite'):
            V.trait_ridge(features, bad, split=SPLIT, exclusions=EXCLUSIONS)
    with pytest.raises(ValueError, match='Infinite features'):
        V.trait_ridge(features.assign(x=np.inf), values, split=SPLIT, exclusions=EXCLUSIONS)
    with pytest.raises(ValueError, match='unique frozen'):
        V.trait_ridge(pd.concat([features,features[['x']]],axis=1), values, split=SPLIT, exclusions=EXCLUSIONS)


def test_unsupported_features_keep_every_test_row_as_an_abstention():
    features, values = inputs()
    batch = V.trait_ridge(features*0+np.nan, values, split=SPLIT, exclusions=EXCLUSIONS)
    assert batch.rows.entity.tolist()==list(SPLIT.entities('test'))
    assert batch.rows.abstained.all() and batch.rows.prediction.isna().all()
    assert batch.model_state['coefficients'] is None
