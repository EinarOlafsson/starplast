"""Bounded train-fit numeric operators for typed precomputation jobs.

This intermediate artifact contains permitted feature coordinates and their
training distributions, never labels, biological predictions or confidence.
Training ranks follow the existing native convention; held-out values use only
the frozen training ECDF. Joint graph or representation construction is outside
this partition. Source snapshots and exclusions remain exact cache dependencies.
"""
from __future__ import annotations

from dataclasses import asdict,dataclass
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from . import artifacts as A,precompute_jobs as P,strategies as S,splits as SP
from .splits import ExclusionManifest,SplitManifest

PHASE='intermediate_numeric_features'
TRANSFORM='training_average_rank_percentile_centered; queries_frozen_training_ecdf'
_ATTENTION=re.compile(r'^(?:n_)?(?:publications?|full_?texts?|citations?|papers?)(?:_|$)',re.IGNORECASE)


@dataclass(frozen=True)
class FeatureLimits:
    """Explicit bounds checked before reading/parsing or allocating numeric matrices."""

    max_rows: int
    max_columns: int
    max_cells: int
    max_input_bytes: int

    def __post_init__(self):
        if any(type(value) is not int or value<=0 for value in asdict(self).values()):
            raise ValueError('Numeric feature limits must be positive integers')

    def check_shape(self,rows,columns):
        """Reject populations exceeding any declared matrix bound."""
        if rows>self.max_rows or columns>self.max_columns or rows*columns>self.max_cells:
            raise P.BudgetExceeded('Numeric feature shape exceeds the declared partition bounds')


def settings(limits):
    """Return the immutable recipe metadata required by this intermediate builder."""
    if not isinstance(limits,FeatureLimits):raise TypeError('Typed numeric feature bounds required')
    return {'phase':PHASE,'feature_transform':TRANSFORM,'missing_feature_value':0.,
        'feature_limits':asdict(limits),'inference_outputs':False}


def exclusions_identity(exclusions):
    """Content identity of the exact training-derived source exclusion manifest."""
    if not isinstance(exclusions,ExclusionManifest):raise TypeError('Typed source exclusions required')
    return hashlib.sha256(A.canonical_object(asdict(exclusions)).encode()).hexdigest()


def code_snapshot():
    """The exact local builder source dependency, suitable for Plan snapshot checks."""
    path=Path(__file__).resolve()
    return P.Snapshot(A.Dependency('code','numeric_feature_builder',hashlib.sha256(path.read_bytes()).hexdigest()),str(path))


def _guard_snapshots():
    return tuple(P.Snapshot(A.Dependency('code','numeric_guard_'+name,
        hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()),str(Path(module.__file__).resolve()))
        for name,module in (('strategies',S),('splits',SP),('artifacts',A),('precompute_jobs',P)))


def snapshot_payload(features, *, organism):
    """Encode an already selected numeric table for caller-owned immutable freezing.

    This returns data only; it writes no file and selects no features or sources.
    Missing values become explicit null and all finite floats retain binary64
    round-trip precision. The builder later verifies both this data and its bytes.
    """
    _validate_frame(features)
    numeric=features.to_numpy(dtype=float,na_value=np.nan)
    if np.isinf(numeric).any():raise ValueError('Infinite features require an explicit preprocessing decision')
    return {'schema_version':1,'organism':organism,'entity_order':list(features.index),
        'columns':list(features.columns),'values':[[float(value) if np.isfinite(value) else None for value in row] for row in numeric]}


def _validate_frame(features):
    if (not isinstance(features,pd.DataFrame) or features.index.has_duplicates or features.columns.has_duplicates
            or not all(isinstance(value,str) and value and value==value.strip() for value in (*features.index,*features.columns))):
        raise ValueError('Numeric inputs require explicit ordered unique string entities and columns')
    if any(pd.api.types.is_complex_dtype(features[column]) for column in features):
        raise ValueError('Complex features require an explicit real-valued preprocessing decision')
    if any(not pd.api.types.is_numeric_dtype(features[column]) and not pd.api.types.is_bool_dtype(features[column]) for column in features):
        raise ValueError('Permitted feature inputs must already have numeric dtypes')


class NumericFeatureBuilder:
    """A cooperative precompute_jobs builder bound to one explicit feature snapshot.

    Supply a Snapshot containing ``snapshot_payload`` JSON, the corresponding
    dataframe, a frozen split/exclusion manifest and allocation bounds. Include
    ``snapshots`` in the Plan. Its reusable_base spec must pin table/code/split/
    exclusions, complete entity order, train-only fitting and ``settings(limits)``.
    The spec's declared output describes the downstream strategy, not a new
    inference kind produced by this intermediate artifact.
    """

    def __init__(self, features, snapshot, *, split, exclusions, limits):
        _validate_frame(features)
        if not isinstance(snapshot,P.Snapshot) or snapshot.dependency.kind!='table':
            raise TypeError('An explicit typed feature-table snapshot is required')
        if not isinstance(split,SplitManifest) or not isinstance(exclusions,ExclusionManifest):
            raise TypeError('Frozen split and source exclusions are required')
        if not isinstance(limits,FeatureLimits):raise TypeError('Explicit typed feature bounds required')
        limits.check_shape(len(features),len(features.columns))
        self.features=features.copy(deep=True)
        self.snapshot=snapshot;self.split=split;self.exclusions=exclusions;self.limits=limits
        self.code=code_snapshot()
        self.guard_code=_guard_snapshots()

    def availability(self):
        """Training-only operator availability, for declaring the immutable job status.

        An unavailable spec needs these gaps before it enters the typed job graph;
        the runner can then retain this partition and block its dependents.
        """
        self._guard_binding()
        self.split.guard_fit('feature_selection',self.split.entities('train'))
        present=self.features.loc[list(self.split.entities('train'))].notna().any().any()
        return {'status':'ok' if present else 'unavailable','gaps':[] if present else
            ['No training-observed numeric features; every input entity retained with an empty operator']}

    def _guard_binding(self):
        ids=tuple(assignment.entity for assignment in self.split.assignments)
        if tuple(self.features.index)!=ids:raise ValueError('Numeric feature population/order differs from the frozen split')
        if self.exclusions.fit_entities!=self.split.entities('train'):
            raise ValueError('Numeric feature fitting requires the complete ordered training cohort')
        if self.exclusions.split_identity!=self.split.identity or self.exclusions.benchmark_id!=self.split.benchmark_id:
            raise ValueError('Numeric feature exclusion/split identities differ')
        self.exclusions.guard_inputs(columns=self.features.columns,benchmark_id=self.split.benchmark_id)

    @property
    def snapshots(self):
        return (self.snapshot,self.code,*self.guard_code)

    def __call__(self, context):
        """Fit numeric ranks on training rows only and preserve every declared query row."""
        if not isinstance(context,P.Context):raise TypeError('Typed precomputation builder context required')
        context.checkpoint()
        spec=context.expected
        if (spec.role!='reusable_base' or spec.fit_role!='train' or spec.evaluation_partition!='none'
                or spec.confidence_kind!='none' or spec.status not in {'ok','unavailable'}):
            raise ValueError('Numeric preprocessing requires a train-fit reusable_base with no inference confidence')
        spec.validate_split(self.split)
        ids=tuple(assignment.entity for assignment in self.split.assignments)
        if spec.entity_order!=ids or tuple(self.features.index)!=ids:
            raise ValueError('Numeric feature population/order differs from the frozen split')
        if spec.fit_entities!=self.split.entities('train') or self.exclusions.fit_entities!=self.split.entities('train'):
            raise ValueError('Numeric feature fitting requires the complete ordered training cohort')
        if (self.exclusions.split_identity!=self.split.identity or self.exclusions.benchmark_id!=self.split.benchmark_id
                or spec.benchmark_id!=self.split.benchmark_id or spec.target not in self.exclusions.targets):
            raise ValueError('Numeric feature target/exclusion/split identities differ')
        if any(snapshot.dependency not in spec.dependencies for snapshot in self.snapshots):
            raise ValueError('Numeric source/code snapshot differs from declared dependencies')
        if [dependency.sha256 for dependency in spec.dependencies if dependency.kind=='exclusions']!=[exclusions_identity(self.exclusions)]:
            raise ValueError('Numeric source closure differs from the declared exclusion dependency')
        if json.loads(spec.settings_json)!=settings(self.limits):
            raise ValueError('Numeric operator recipe or bounds differ from declared settings')
        self.exclusions.guard_inputs(columns=self.features.columns,benchmark_id=self.split.benchmark_id)
        if any(S.NEVER_FEATURES.search(column) or _ATTENTION.search(column) for column in self.features):
            raise ValueError('Identifiers, free text and literature attention cannot become numeric features')
        for stage in ('feature_selection','imputation','scaling'):
            self.split.guard_fit(stage,spec.fit_entities)
        self.limits.check_shape(len(self.features),len(self.features.columns))
        if Path(self.snapshot.path).stat().st_size>self.limits.max_input_bytes:
            raise P.BudgetExceeded('Numeric source snapshot exceeds the declared input byte bound')
        for snapshot in self.snapshots:snapshot.verify()
        source=json.loads(Path(self.snapshot.path).read_bytes())
        expected_source=snapshot_payload(self.features,organism=self.split.organism)
        if A.canonical_object(source)!=A.canonical_object(expected_source):
            raise ValueError('Numeric dataframe does not exactly match its pinned source snapshot')
        numeric=pd.DataFrame(self.features.to_numpy(dtype=float,na_value=np.nan),index=self.features.index,columns=self.features.columns)
        train=numeric.loc[list(spec.fit_entities)]
        kept=train.columns[train.notna().any()].tolist()
        withheld=train.columns[~train.notna().any()].tolist()
        availability=self.availability()
        if spec.status!=availability['status'] or (not kept and not spec.gaps):
            raise ValueError('Declared numeric operator availability differs from training-only evidence')
        matrix=np.zeros((len(numeric),len(kept)),dtype=float)
        available=np.zeros_like(matrix,dtype=bool)
        distributions={}
        for j,column in enumerate(kept):
            context.checkpoint()
            ordered=np.sort(train[column].dropna().to_numpy(dtype=float))
            distributions[column]=ordered.tolist()
            values=numeric[column].to_numpy(dtype=float)
            known=np.isfinite(values);available[:,j]=known
            left=np.searchsorted(ordered,values[known],side='left')
            right=np.searchsorted(ordered,values[known],side='right')
            ranks=np.where(right>left,(left+right+1)/2,right)/len(ordered)
            matrix[known,j]=ranks-.5
        context.checkpoint()
        for snapshot in self.snapshots:snapshot.verify()
        return {'operator.json':{'schema_version':1,'phase':PHASE,'status':'available' if kept else 'unavailable',
            'entity_order':list(ids),'columns':kept,'input_columns':list(numeric.columns),
            'withheld_all_missing_training_columns':withheld,'matrix':matrix.tolist(),'observed':available.tolist(),
            'training_distributions':distributions,'fit_entities':list(spec.fit_entities),'fit_role':'train',
            'split_identity':self.split.identity,'benchmark_id':self.split.benchmark_id,
            'source_sha256':self.snapshot.dependency.sha256,'exclusions_sha256':exclusions_identity(self.exclusions),
            'source_exclusion_columns':list(self.exclusions.columns),'source_exclusion_layers':list(self.exclusions.layers),
            'feature_transform':TRANSFORM,'missing_feature_value':0.,
            'gaps':availability['gaps'],
            'inference_outputs':False,'biological_accuracy':None,'calibrated_confidence':None}}
