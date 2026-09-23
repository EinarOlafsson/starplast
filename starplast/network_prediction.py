"""Family-held-out diffusion over explicitly sourced network layers.

This method is transductive: all graph nodes and fixed edges are visible, but test
labels are never seeds or tuning targets. Layer weights and temperature use two
different inner validation sets. Unsupported nodes explicitly abstain.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib

import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.special import softmax
from sklearn.metrics import log_loss

from .prediction import (PredictionResult, TaskSpec, _evaluate, _splits, group_ids,
                         known_labels)


def run_network(nodes, spec, graph, layers, edge_sources, log=print):
    """Evaluate weighted propagation with explicit edge provenance and aligned IDs.

    ``graph`` must contain ordered ``gene_ids`` and the usual ``layer__a/b/w``
    arrays. ``edge_sources`` maps every layer to the columns used to build it;
    registered target-derived layers are excluded even if omitted from this map.
    User-supplied edge provenance must describe all label-derived inputs.
    """
    from .methods import Propagator, layer_matrix
    from .search import excluded_detail, excluded_edges
    from .datasets import derived_dependents
    if spec.kind != 'classification':
        raise ValueError('network propagation currently supports classification')
    spec = replace(spec, method="network")
    nodes = nodes.reset_index(drop=True)
    if nodes.gene_id.isna().any() or nodes.gene_id.duplicated().any():
        raise ValueError("one unique nonmissing gene ID per row is required")
    if "gene_ids" not in graph or not np.array_equal(np.asarray(graph['gene_ids']).astype(str), nodes.gene_id.astype(str)):
        raise ValueError("network ordered gene IDs do not match the input table")
    missing = set(layers)-set(edge_sources)
    if missing:
        raise ValueError(f"edge source provenance is required for: {sorted(missing)}")
    # Static closure only: no test labels participate in layer availability.
    hidden = nodes.copy()
    hidden[spec.target] = None
    banned_columns = set(excluded_detail(hidden, spec.target)) | {spec.target} | set(spec.exclude)
    banned_columns |= set(derived_dependents(banned_columns))
    banned_layers = excluded_edges(spec.target, scope="target_family")
    excluded = {name: "target-derived graph or source" for name in layers
                if name in banned_layers or set(edge_sources[name]) & banned_columns}
    layers = [name for name in layers if name not in excluded]
    if not layers:
        raise ValueError("all network layers derive from the held-out target")
    operators = {}
    for name in layers:
        for key in ('a','b','w'):
            if f'{name}__{key}' not in graph:
                raise ValueError(f'missing network array: {name}__{key}')
        weights = np.asarray(graph[f'{name}__w'])
        a,b = (np.asarray(graph[f'{name}__{part}']) for part in ('a','b'))
        if not (a.ndim == b.ndim == weights.ndim == 1 and len(a)==len(b)==len(weights)):
            raise ValueError('network endpoint and weight arrays must have matching lengths')
        if any(endpoint.dtype.kind not in 'iu' or (endpoint<0).any() or (endpoint>=len(nodes)).any()
               for endpoint in (a,b)):
            raise ValueError('network endpoints must be integer indices into the gene table')
        if not np.isfinite(weights).all() or (weights < 0).any():
            raise ValueError("diffusion requires finite non-negative weights")
        operators[name] = layer_matrix(name,len(nodes),graph)
    known = known_labels(nodes[spec.target])
    classes = sorted(set(nodes.loc[known,spec.target].map(str)))
    if len(classes)<2:
        raise ValueError("at least two observed classes required")
    lookup = {label:i for i,label in enumerate(classes)}
    y = np.array([lookup.get(str(v),-1) for v in nodes[spec.target]])
    groups = group_ids(nodes,spec.group_column)
    indices = np.flatnonzero(known)
    if len(indices)<10:
        raise ValueError('at least ten labelled genes are required for evaluation')
    n, c = len(nodes), len(classes)

    def fit(train):
        """Separate seeding, layer-selection and calibration groups."""
        base, calibration = _splits(train,y,groups,spec,n_splits=5)[0]
        seeds, validation = _splits(base,y,groups,spec,n_splits=4)[0]
        fields=[]
        for name in layers:
            model=Propagator(operators[name]).fit(seeds[:,None],y[seeds])
            scores=np.zeros((n,c))
            scores[:,model.classes_]=model.scores_.T
            total=scores.sum(axis=1,keepdims=True)
            scores=np.divide(scores,total,out=np.zeros_like(scores),where=total>0)
            fields.append(scores)
        fields=np.stack(fields,axis=0)
        # Optimize a convex mixture against inner validation labels, never outer test labels.
        def mix(weights):
            score=np.einsum('l,lnc->nc',weights,fields)
            total=score.sum(axis=1,keepdims=True)
            return np.divide(score,total,out=np.full_like(score,1/c),where=total>0)
        def loss(weights):
            return log_loss(y[validation],mix(weights)[validation],labels=np.arange(c))+.01*np.sum(weights**2)
        initial=np.full(len(layers),1/len(layers))
        if len(layers)>1:
            opt=minimize(loss,initial,method='SLSQP',bounds=[(0,1)]*len(layers),
                         constraints={'type':'eq','fun':lambda w: w.sum()-1})
            if not opt.success:
                raise RuntimeError(f"network weight optimization failed: {opt.message}")
            weights=np.clip(opt.x,0,1);weights/=weights.sum()
        else:
            weights=initial
        probabilities=mix(weights)
        support=np.einsum('l,lnc->nc',weights,fields).sum(axis=1)>1e-12
        temperature=1.
        status='not_calibrated'
        eligible=calibration[support[calibration]]
        if spec.calibrate and len(eligible)>=max(20,2*c) and len(np.unique(y[eligible]))>1:
            logs=np.log(np.clip(probabilities,1e-12,1))
            def loss_temperature(value):
                return log_loss(y[eligible],softmax(logs[eligible]/np.exp(value),axis=1),labels=np.arange(c))
            opt=minimize_scalar(loss_temperature,bounds=(-2,3),method='bounded')
            temperature=float(np.exp(opt.x))
            probabilities=softmax(logs/temperature,axis=1)
            status='held_out_temperature_scaling'
        elif spec.calibrate:
            status='insufficient_calibration_labels'
        report={'training_ids':nodes.gene_id.iloc[seeds].tolist(),
                'weight_validation_ids':nodes.gene_id.iloc[validation].tolist(),
                'calibration_ids':nodes.gene_id.iloc[calibration].tolist(),
                'layer_weights':dict(zip(layers,weights.tolist())),
                'temperature':temperature,'calibration':status,'excluded_layers':excluded}
        return probabilities,support,report

    scores=np.zeros((n,c));supported=np.zeros(n,dtype=bool);fold_ids=np.full(n,-1)
    reports=[];statuses=np.full(n,'not_calibrated',dtype=object)
    for fold,(train,test) in enumerate(_splits(indices,y,groups,spec)):
        log(f'Network fold {fold+1}/{spec.folds}')
        probabilities,support,report=fit(train)
        scores[test]=probabilities[test];supported[test]=support[test];fold_ids[test]=fold
        statuses[test]=report['calibration']
        report.update(fold=fold,test_ids=nodes.gene_id.iloc[test].tolist())
        reports.append(report)
    probabilities,support,final=fit(indices)
    scores[~known]=probabilities[~known];supported[~known]=support[~known]
    statuses[~known]=final['calibration']
    supported &= scores.max(axis=1)>=spec.min_probability
    predicted=scores.argmax(axis=1)
    metrics,per_class=_evaluate(y,predicted,scores,known,supported,classes,'classification')
    table=pd.DataFrame({'gene_id':nodes.gene_id,'group':groups,'fold':fold_ids,
                        'role':np.where(known,'held_out_evaluation','unlabelled_candidate'),
                        'truth':nodes[spec.target],'prediction':[classes[i] for i in predicted],
                        'supported':supported,'score':scores.max(axis=1),
                        'calibration_status':statuses,'evidence_status':'model_prediction',
                        'abstention_reason':np.where(supported,'','no_seed_support_or_low_probability')})
    table.loc[~supported,'prediction']=None
    for i,label in enumerate(classes):
        table[f'probability::{label}']=scores[:,i]
    fingerprints={name:hashlib.sha256(b''.join(np.asarray(graph[f'{name}__{part}']).tobytes()
                                             for part in ('a','b','w'))).hexdigest() for name in layers}
    provenance={'evaluation_mode':'transductive_group_holdout','ordered_gene_ids':nodes.gene_id.tolist(),
                'edge_sources':edge_sources,'edge_sha256':fingerprints,'final_fit':final,
                'executed_method':'weighted_network_diffusion','classes':classes,
                'limitations':['All fixed graph nodes are visible. Performance is not an inductive new-node estimate.',
                               'Edge provenance must exclude all target-derived information.']}
    return PredictionResult(spec,table,metrics,per_class,reports,provenance)
