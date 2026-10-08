"""Frozen categorical ortholog transfer with receiver-training-only mapping.

Donor aggregation and orthogroup projection precede this adapter and require
their own source/mapping identities. Missing mappings and unsupported donor
values remain abstentions, rather than shrinking the evaluation denominator.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import re

import numpy as np
import pandas as pd

from . import organisms as O, strategy_catalog as C
from .label_records import LabelBatch
from .provenance import EVIDENCE_GRADES
from .splits import ExclusionManifest, SplitManifest


@dataclass(frozen=True)
class TransferSource:
    """Explicit donor lineage, mapping identity and semantic/dependence review."""

    organism: str
    column: str
    table_sha256: str
    mapping_sha256: str
    evidence_grade: str
    compatibility_note: str
    receiver_truth_dependency: str = 'unknown'

    def __post_init__(self):
        if self.organism not in O.SPACES:
            raise ValueError('Donor needs an implemented organism space')
        if not all(isinstance(v, str) and v.strip() for v in (self.column, self.evidence_grade, self.compatibility_note)):
            raise ValueError('Source column, evidence grade and semantic review must be explicit')
        if self.evidence_grade not in EVIDENCE_GRADES:
            raise ValueError('Use a registered donor evidence grade')
        for digest in (self.table_sha256, self.mapping_sha256):
            if not isinstance(digest, str) or not re.fullmatch(r'[a-f0-9]{64}', digest):
                raise ValueError('Donor table and mapping need SHA-256 identities')
        if self.receiver_truth_dependency not in {'unknown', 'independent', 'derived_receiver'}:
            raise ValueError('Declare donor dependence on receiver truth')


def ortholog_transfer(mapped_source, training_labels, *, split, exclusions, source, numeric_source):
    """Retain native source-to-label calls for every frozen outer-test entity.

    ``mapped_source`` contains donor values already projected onto the ordered
    split population through the native orthogroup aggregator. The native
    transfer function fits categorical modes or numeric quantile bins on
    receiver training genes only. It produces no native class-score matrix or
    probability; support and calibrated confidence remain unavailable.
    """
    if not isinstance(split, SplitManifest) or not isinstance(exclusions, ExclusionManifest):
        raise TypeError('Frozen split and training-derived exclusions are required')
    if exclusions.split_identity != split.identity or exclusions.benchmark_id != split.benchmark_id:
        raise ValueError('Source exclusions and split identities differ')
    if exclusions.fit_entities != split.entities('train'):
        raise ValueError('Exclusions must use the complete ordered training cohort')
    if not isinstance(source, TransferSource) or source.organism == split.organism:
        raise ValueError('An explicitly different donor species is required')
    if source.receiver_truth_dependency == 'derived_receiver':
        raise ValueError('Donor values derived from receiver truth contaminate evaluation')
    ids = tuple(row.entity for row in split.assignments)
    if not isinstance(mapped_source, pd.Series) or tuple(mapped_source.index) != ids:
        raise ValueError('Mapped donor population/order differs from the frozen split')
    if not isinstance(training_labels, pd.Series) or tuple(training_labels.index) != split.entities('train'):
        raise ValueError('Provide the complete ordered receiver training-label cohort only')
    if not all(isinstance(v, str) and v for v in training_labels):
        raise ValueError('Receiver training labels must be observed categorical values')
    if type(numeric_source) is not bool:
        raise ValueError('Declare numeric versus categorical donor semantics')
    if numeric_source:
        values = pd.to_numeric(mapped_source, errors='raise').astype(float)
        if np.isinf(values.to_numpy()).any():
            raise ValueError('Infinite donor values require an explicit source decision')
    else:
        values = mapped_source.astype(object)
        if not all(isinstance(v, str) and v for v in values.dropna()):
            raise ValueError('Known categorical donor values must be nonempty strings')
    # Receiver column bans cannot be applied by name to a different species.
    # Donor dependence is recorded separately; known receiver-derived truth is refused.
    for stage in ('feature_selection', 'scaling', 'model_fit'):
        split.guard_fit(stage, training_labels.index)
    train = values.loc[list(training_labels.index)]
    mapping = C._transfer_map(train, training_labels, numeric_source, False)
    native = mapping(values).where(values.notna())
    test = list(split.entities('test'))
    predictions = native.loc[test].astype(object).where(native.loc[test].notna(), None).tolist()
    rows = pd.DataFrame({'entity': test, 'prediction': pd.Series(predictions, dtype=object),
        'support': [None] * len(test), 'abstained': [not isinstance(v, str) for v in predictions],
        'calibrated_confidence': [None] * len(test),
        'donor_mapped': values.loc[test].notna().to_numpy()})
    model = {'strategy': 'ortholog_transfer', 'split_identity': split.identity,
        'benchmark_id': split.benchmark_id, 'fit_role': 'train', 'fit_entities': list(training_labels.index),
        'source': asdict(source), 'numeric_source': numeric_source,
        'training_mapped_source': train.astype(object).where(train.notna(), None).tolist(),
        'training_labels': training_labels.tolist(), 'training_class_counts': training_labels.value_counts().to_dict(),
        'training_mapped_entities': train.index[train.notna()].tolist(),
        'source_exclusion_columns': list(exclusions.columns), 'source_exclusion_layers': list(exclusions.layers),
        'native_class_scores': 'unavailable', 'native_support': 'unavailable'}
    gaps = ['Native ortholog transfer produces no class scores or calibrated confidence',
        'Biological accuracy requires admitted receiver truth and compatible donor context',
        'Unmapped and unsupported donor values are retained as explicit abstentions']
    if source.receiver_truth_dependency == 'unknown':
        gaps.append('Donor independence from receiver truth is unresolved')
    if not train.notna().any():
        gaps.append('No training donor mapping; every test entity retained as an abstention')
    return LabelBatch(rows, pd.DataFrame(index=test), model, tuple(gaps))
