"""Frozen native absolute-residual conformal intervals from train-only numeric models."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import artifacts as A, scorecard as SC, strategy_learning as N
from .splits import SplitManifest
from .value_records import ValueBatch


def conformal_values(calibration_predictions, test_predictions, calibration_values, *, split, base_model, alpha=.1):
    """Calibrate only on the complete separate calibration cohort; retain all test rows.

    The verified base model must fit only training values. Missing calibration
    predictions do not shrink the calibration population: intervals become
    explicitly unavailable. Infinite native quantiles have typed unbounded status,
    never nonfinite JSON. No outer-test truth enters this adapter.
    """
    if not isinstance(split, SplitManifest) or not isinstance(base_model, A.Artifact):
        raise TypeError('Frozen split and verified fitted-model artifact are required')
    spec = base_model.spec
    spec.validate_split(split)
    if (spec.role != 'fitted_model' or spec.strategy != 'trait_regression' or spec.task != SC.T_VALUES
            or spec.output.kind != 'numeric_estimates'):
        raise ValueError('Use a native numeric regression fitted-model artifact')
    if spec.fit_role != 'train' or spec.fit_entities != split.entities('train') or spec.benchmark_id != split.benchmark_id:
        raise ValueError('Base numeric model must fit only the complete ordered training cohort')
    base_model.verify_contents()
    state = base_model.payloads.get('model_state.json', {})
    if state.get('fit_entities') != list(spec.fit_entities):
        raise ValueError('Base model payload and fitting population differ')
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not 0 < alpha < 1:
        raise ValueError('Declare a finite conformal alpha strictly between zero and one')
    for values, role in ((calibration_predictions,'calibration'),(test_predictions,'test'),(calibration_values,'calibration')):
        if not isinstance(values,pd.Series) or tuple(values.index) != split.entities(role):
            raise ValueError('Provide the complete ordered '+role+' cohort only')
        if np.isinf(values.to_numpy(dtype=float)).any():
            raise ValueError('Infinite predictions/values require an explicit source decision')
    truth = calibration_values.to_numpy(dtype=float)
    if not np.isfinite(truth).all():
        raise ValueError('Calibration values must be explicitly observed finite measurements')
    split.guard_fit('calibration',calibration_values.index)
    calibration = calibration_predictions.to_numpy(dtype=float)
    prediction = test_predictions.to_numpy(dtype=float)
    with np.errstate(over='ignore'):
        residual = np.abs(calibration-truth)
    if np.isinf(residual).any():
        raise ValueError('Finite calibration residuals overflowed the declared numeric representation')
    supported = np.isfinite(residual).all()
    half = N._conformal_quantile(residual,alpha) if supported else np.nan
    threshold = {'status':'finite' if np.isfinite(half) else 'unbounded' if np.isposinf(half) else 'unavailable',
                 'value':float(half) if np.isfinite(half) else None}
    status = [threshold['status'] if np.isfinite(value) else 'unavailable' for value in prediction]
    rows = pd.DataFrame({'entity':list(split.entities('test')),
        'prediction':[float(value) if np.isfinite(value) else None for value in prediction],
        'lower':[float(value-half) if s=='finite' else None for value,s in zip(prediction,status)],
        'upper':[float(value+half) if s=='finite' else None for value,s in zip(prediction,status)],
        'interval_status':status,'abstained':~np.isfinite(prediction),
        'calibrated_confidence':[None]*len(prediction)})
    if not np.isfinite(rows.loc[rows.interval_status.eq('finite'),['lower','upper']].to_numpy(dtype=float)).all():
        raise ValueError('Native finite interval bounds overflowed the declared numeric representation')
    calibration_state = {'strategy':'conformal_values','split_identity':split.identity,'benchmark_id':split.benchmark_id,
        'base_model_identity':base_model.identity,'fit_entities':list(spec.fit_entities),'fit_role':'train',
        'calibration_entities':calibration_values.index.tolist(),'calibration_values':truth.tolist(),
        'calibration_predictions':[float(value) if np.isfinite(value) else None for value in calibration],
        'absolute_residuals':[float(value) if np.isfinite(value) else None for value in residual],
        'alpha':alpha,'half_width':threshold,'missing_calibration_predictions':int((~np.isfinite(calibration)).sum())}
    gaps = ['Empirical interval coverage requires independent source/context admission and exchangeability review',
            'Nominal coverage is not a measured per-gene probability or a deployment guarantee',
            'One fixed numeric base/cohort; other variants and full outer coverage remain open']
    if not supported:
        gaps.append('Calibration predictions incomplete; no reduced-population quantile or interval fabricated')
    if threshold['status']=='unbounded':
        gaps.append('Insufficient calibration quantile support; unbounded intervals have no finite inference capacity')
    return ValueBatch(rows,calibration_state,tuple(gaps))
