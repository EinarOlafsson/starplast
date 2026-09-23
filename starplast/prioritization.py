"""Transparent candidate explanations and experiment-budget prioritization.

Predictions remain hypotheses. Expected value uses explicit user benefits/costs;
uncertainty-driven exploration is reported separately from success probability.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def explain_candidates(result, nodes, limit=50):
    """Explain supported unknown-gene calls with observed sources and training neighbours.

    Linear contributions describe the fitted model's logits before temperature
    scaling. For nonlinear models, neighbours provide context rather than causal
    feature attributions. Held-out evaluation rows are never explained using the
    final refit, which has seen their labels.
    """
    from .datasets import provenance
    from sklearn.neighbors import NearestNeighbors
    fitted=result.model
    if fitted is None:
        raise ValueError('this result has no fitted feature model for explanations')
    data=nodes.set_index('gene_id',drop=False)
    candidates=result.predictions.query("role == 'unlabelled_candidate' and supported").copy()
    if 'score' in candidates: candidates=candidates.sort_values('score',ascending=False)
    candidates=candidates.head(limit)
    train=data.loc[result.provenance['final_training_ids']]
    representation=fitted._matrix(train)
    neighbors=NearestNeighbors(n_neighbors=min(5,len(train))).fit(representation)
    explained=[]
    for row in candidates.itertuples(index=False):
        gene=data.loc[[row.gene_id]]
        distances,positions=neighbors.kneighbors(fitted._matrix(gene))
        sources={}
        for feature in fitted.transformer.columns_:
            value=gene[feature].iloc[0]
            if pd.isna(value): continue
            origin=provenance(feature)
            key='esm2_protein' if feature.startswith('esm_') else origin.key if origin else 'unregistered'
            sources[key]=sources.get(key,0)+1
        contributions=[]
        if result.spec.method=='linear' and hasattr(fitted.model,'coef_'):
            coefficients=fitted.model.coef_
            if np.ndim(coefficients)==1:
                coefficient=coefficients
            elif len(coefficients)==1:
                positive=fitted.classes[int(fitted.model.classes_[1])]
                coefficient=coefficients[0]*(1 if str(row.prediction)==str(positive) else -1)
            else:
                class_index=fitted.classes.index(str(row.prediction))
                coefficient=coefficients[list(fitted.model.classes_).index(class_index)]
            values=fitted.transformer.transform(gene)[0]*coefficient
            for j in np.argsort(np.abs(values))[-8:][::-1]:
                feature=fitted.transformer.columns_[j]
                if pd.isna(gene[feature].iloc[0]): continue
                source=provenance(feature)
                contributions.append({'feature':feature,'value':float(gene[feature].iloc[0]),
                                      'model_contribution':float(values[j]),
                                      'contribution_scale':'logit' if result.spec.kind=='classification' else 'outcome',
                                      'source':source.key if source else 'esm2_protein' if feature.startswith('esm_') else None,
                                      'source_url':source.url if source else None})
        explained.append({'gene_id':row.gene_id,'prediction':row.prediction,'measured_features_by_source':sources,
                          'training_neighbors':[{'gene_id':str(train.gene_id.iloc[j]),
                                                'label':str(train[result.spec.target].iloc[j]),
                                                'distance':float(distance)}
                                               for j,distance in zip(positions[0],distances[0])],
                          'linear_contributions':contributions,
                          'caution':'Model associations suggest experiments; they do not establish mechanism.'})
    return explained


def prioritize(candidates, budget, probability_column='score', cost_column='cost',
               benefit_column='benefit', group_column=None, exploration_weight=0.):
    """Select a greedy budget-constrained list, at most one candidate per supplied group.

    Utility is ``p*benefit + exploration_weight*4*p*(1-p)`` per cost. This heuristic
    is not an optimal knapsack solver. Probabilities should come from separately
    calibrated predictions; score rankings alone must not be passed as probabilities.
    Missing cost/benefit defaults to one and is recorded in the output.
    """
    if budget<=0 or not 0<=exploration_weight<=1:
        raise ValueError('positive budget and exploration_weight in [0,1] required')
    table=candidates.copy()
    if 'supported' in table: table=table[table.supported.fillna(False)].copy()
    if 'calibration_status' not in table or not table.calibration_status.eq('held_out_temperature_scaling').all():
        raise ValueError('expected-success prioritization requires calibrated probabilities')
    if table.gene_id.duplicated().any(): raise ValueError('duplicate candidate gene IDs')
    probability=pd.to_numeric(table[probability_column],errors='coerce')
    if not probability.between(0,1).all(): raise ValueError('probabilities must be finite and between zero and one')
    for column in [cost_column,benefit_column]:
        defaulted=column not in table
        if defaulted: table[column]=1.
        table[column+'_assumption']='default one' if defaulted else 'user supplied'
        table[column]=pd.to_numeric(table[column],errors='raise')
    if not np.isfinite(table[[cost_column,benefit_column]].to_numpy()).all():
        raise ValueError('cost and benefit must be finite')
    if (table[cost_column]<=0).any() or (table[benefit_column]<0).any():
        raise ValueError('cost must be positive and benefit non-negative')
    table['expected_benefit']=probability*table[benefit_column]
    table['exploration_value']=exploration_weight*4*probability*(1-probability)
    table['utility_per_cost']=(table.expected_benefit+table.exploration_value)/table[cost_column]
    table=table.sort_values(['utility_per_cost','gene_id'],ascending=[False,True])
    selected=[];spent=0.;groups=set()
    for index,row in table.iterrows():
        group=row[group_column] if group_column else None
        if group_column and pd.isna(group): group='gene:'+str(row.gene_id)
        if group_column and group in groups: continue
        if spent+row[cost_column] > budget+1e-12: continue
        selected.append(index);spent+=row[cost_column]
        if group_column: groups.add(group)
    result=table.loc[selected].copy()
    result['cumulative_cost']=result[cost_column].cumsum()
    result['selection_method']='greedy expected benefit plus stated exploration bonus'
    return result


def compare_screen(nodes, screen, value_column, gene_column='gene_id'):
    """Join a measured screen to gene evidence without dropping unmapped identifiers.

    Duplicate screen genes are refused so guide-level replicates are not silently
    averaged. Summarize guides in spaCR first, then import the gene-level effects.
    """
    if gene_column not in screen or value_column not in screen:
        raise ValueError('screen needs a gene identifier and selected effect column')
    if screen[gene_column].duplicated().any():
        raise ValueError('summarize duplicate gene rows before comparing this screen')
    from .identity import build_index, norm
    from .paths import cache_file
    species='Pf' if nodes.gene_id.astype(str).str.startswith('PF3D7').any() else 'Tg'
    identity=build_index(nodes.gene_id,cache_file('plasmodb_identity.tsv' if species=='Pf' else 'toxodb_identity.tsv'),
                         log=lambda *a:None)
    result=screen.copy()
    canonical=set(nodes.gene_id)
    result['input_gene_id']=result[gene_column].astype(str)
    def resolve(value):
        if value in canonical:return value
        hit=identity.lookup.get(norm(value))
        return hit[0] if hit else None
    result['gene_id']=result.input_gene_id.map(resolve)
    mapped=result.gene_id.dropna()
    if mapped.duplicated().any():
        raise ValueError('multiple input identifiers resolve to one gene; resolve replicates explicitly')
    result['mapping_status']=np.where(result.gene_id.notna(),'mapped','unresolved')
    result['screen_effect']=pd.to_numeric(result[value_column],errors='coerce')
    result.loc[~np.isfinite(result.screen_effect),'screen_effect']=np.nan
    result['measurement_status']=np.where(result.screen_effect.notna(),'measured','missing_or_non_numeric')
    return result.merge(nodes,on='gene_id',how='left',validate='many_to_one',suffixes=('_screen',''))
