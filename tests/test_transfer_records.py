"""Native ortholog-map parity without hidden receiver truth or denominator loss."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from starplast import organisms as O, strategy_catalog as C, transfer_records as T
from starplast.splits import Assignment, ExclusionManifest, SplitManifest

IDS = tuple(f'TGME49_{i:06d}' for i in range(100001, 100011))
ROLES = ['train'] * 4 + ['tune', 'calibration'] + ['test'] * 4
SPLIT = SplitManifest(O.TOXOPLASMA, 'transfer-fixture', 'entity',
    tuple(Assignment(e, e, role) for e, role in zip(IDS, ROLES)), 17)
EXCLUSIONS = ExclusionManifest('transfer-fixture', ('target',), ('target',), (), IDS[:4], SPLIT.identity)
SOURCE = T.TransferSource(O.FALCIPARUM, 'target', 'a' * 64, 'b' * 64, 'prediction',
    'Synthetic source; biological compatibility not established')
LABELS = pd.Series(['A', 'B', 'A', 'B'], index=IDS[:4])


@pytest.mark.parametrize('numeric', [False, True])
def test_native_transfer_parity_and_complete_test_population(numeric):
    values = ([1., 2., 3., 4., 10., None, 1., 3., 8., None] if numeric else
              ['x', 'y', 'x', 'y', 'x', None, 'x', 'y', 'unseen', None])
    mapped = pd.Series(values, index=IDS)
    batch = T.ortholog_transfer(mapped, LABELS, split=SPLIT, exclusions=EXCLUSIONS,
        source=SOURCE, numeric_source=numeric)
    visible = pd.Series([None] * len(IDS), index=IDS)
    visible.loc[LABELS.index] = LABELS
    native = C._transfer_map(mapped, visible, numeric, False)(mapped).where(mapped.notna()).loc[list(IDS[6:])]
    assert batch.rows.prediction.tolist() == native.astype(object).where(native.notna(), None).tolist()
    assert batch.rows.entity.tolist() == list(IDS[6:])
    assert batch.rows.iloc[-1].abstained
    assert batch.class_scores.shape == (4, 0)
    assert batch.rows.support.isna().all() and batch.rows.calibrated_confidence.isna().all()


def test_hidden_donor_distribution_cannot_change_training_map():
    mapped = pd.Series(np.arange(10, dtype=float), index=IDS)
    first = T.ortholog_transfer(mapped, LABELS, split=SPLIT, exclusions=EXCLUSIONS, source=SOURCE, numeric_source=True)
    mapped.loc[list(IDS[4:])] = 1e8
    second = T.ortholog_transfer(mapped, LABELS, split=SPLIT, exclusions=EXCLUSIONS, source=SOURCE, numeric_source=True)
    assert first.model_state == second.model_state


def test_no_mapped_training_genes_abstains_for_every_test_gene():
    mapped = pd.Series([None] * 4 + ['x'] * 6, index=IDS)
    batch = T.ortholog_transfer(mapped, LABELS, split=SPLIT, exclusions=EXCLUSIONS, source=SOURCE, numeric_source=False)
    assert len(batch.rows) == 4 and batch.rows.abstained.all()
    assert any('No training donor' in gap for gap in batch.gaps)


def test_receiver_truth_leaks_and_unknown_source_address_are_refused():
    mapped = pd.Series(['x'] * 10, index=IDS)
    with pytest.raises(ValueError, match='contaminate'):
        T.ortholog_transfer(mapped, LABELS, split=SPLIT, exclusions=EXCLUSIONS,
            source=replace(SOURCE, receiver_truth_dependency='derived_receiver'), numeric_source=False)
    with pytest.raises(ValueError, match='training-label cohort only'):
        T.ortholog_transfer(mapped, pd.Series(['A'] * 10, index=IDS), split=SPLIT, exclusions=EXCLUSIONS,
            source=SOURCE, numeric_source=False)
    with pytest.raises(ValueError, match='population/order'):
        T.ortholog_transfer(mapped.iloc[::-1], LABELS, split=SPLIT, exclusions=EXCLUSIONS,
            source=SOURCE, numeric_source=False)
    with pytest.raises(ValueError, match='SHA-256'):
        replace(SOURCE, mapping_sha256='unresolved')


def test_numeric_source_without_observed_training_mapping_is_explicit_abstention():
    mapped = pd.Series([None] * 4 + [1.] * 6, index=IDS, dtype=float)
    batch = T.ortholog_transfer(mapped, LABELS, split=SPLIT, exclusions=EXCLUSIONS, source=SOURCE, numeric_source=True)
    assert batch.rows.abstained.all()
