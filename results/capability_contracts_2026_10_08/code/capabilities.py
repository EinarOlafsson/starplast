"""Declare strategy meanings and resolve typed questions without running inference.

Applicability is a checked adapter plan, not a performance claim. Catalogue
self-test tasks can differ from deployment outputs (a ranking can rank pairs or
classes; an ablation returns evidence effects). Missing truth and untested
technique roles stay visible. Scores are never converted to gene probabilities.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

from . import organisms as O, scorecard as SC, strategies as S, techniques as T
from .query import BiologicalContext, Query

SCHEMA_VERSION = 1
KINDS = frozenset({'label_calls', 'label_sets', 'numeric_estimates', 'numeric_intervals',
                   'gene_ranking', 'pair_ranking', 'class_recovery', 'map_modules',
                   'gene_set', 'feature_enrichment', 'evidence_effects', 'replication_findings'})
UNITS = frozenset({'gene', 'pair', 'class', 'module', 'feature', 'evidence_family', 'finding'})
ROLES = frozenset({'direct_prediction', 'pipeline_ablation', 'known_truth_null_control'})


@dataclass(frozen=True)
class Output:
    """One allowed output meaning; score semantics are not calibration claims."""

    kind: str
    unit: str
    semantics: str

    def __post_init__(self):
        if self.kind not in KINDS or self.unit not in UNITS or not self.semantics.strip():
            raise ValueError('Output needs an explicit kind, entity unit and score meaning')


@dataclass(frozen=True)
class TechniqueValidation:
    """Required component test, distinct from measured composite performance."""

    technique: str
    roles: tuple[str, ...]
    test: str
    status: str = 'declared_not_measured'

    def __post_init__(self):
        if self.technique not in T.TECHNIQUES or not self.roles or set(self.roles) - ROLES:
            raise ValueError('Unknown technique or validation role')
        if not self.test or self.status != 'declared_not_measured':
            raise ValueError('This declaration cannot assert component accuracy')


@dataclass(frozen=True)
class Capability:
    """Catalogue adapter declaration, including query and fitting dependencies."""

    strategy: str
    query_kinds: tuple[str, ...]
    outputs: tuple[Output, ...]
    required: tuple[str, ...]
    column_parameters: tuple[tuple[str, str], ...]
    target_parameter: str
    minimum_seeds: int
    parameters: tuple[str, ...]
    benchmark_tasks: tuple[str, ...]
    techniques: tuple[TechniqueValidation, ...]
    precomputable: tuple[str, ...]
    cache_dependencies: tuple[str, ...]
    organisms: tuple[str, ...] = (O.TOXOPLASMA, O.FALCIPARUM)

    def validate_output(self, output: Output):
        """Refuse undeclared output kinds, units or reinterpretations of scores."""
        if output not in self.outputs:
            raise ValueError('Adapter cannot change output semantics or manufacture probabilities')

    def to_dict(self):
        """Return a versionable transport record without executable defaults."""
        return asdict(self)


# Each deployment declaration is reviewed independently of its legacy test task.
_CALL = Output('label_calls', 'gene', 'label with method-specific support; not calibrated probability')
_MAP = Output('map_modules', 'module', 'geometry or module membership; not label probability')
_RANK = Output('gene_ranking', 'gene', 'relative method score; not calibrated membership probability')
_PAIR = Output('pair_ranking', 'pair', 'relative relationship support; missing edges are not negatives')
_VALUE = Output('numeric_estimates', 'gene', 'estimated target value in declared target units')
_CLASS = Output('class_recovery', 'class', 'class recovery and coverage; not per-gene probability')
_REPL = Output('replication_findings', 'finding', 'discovery association; independent replication requires disjoint validation')

# key: query kinds, outputs, evidence requirements, typed required columns,
#      primary query target parameter, minimum query seeds.
_SPECS = {
    'holdout_search': ('gene label class', (_MAP, _CALL, _CLASS), 'features', 'target:categorical', 'target', 0),
    'geneset_hunt': ('gene_set', (_MAP, Output('gene_set', 'gene', 'retrieved members; unmeasured membership remains unknown')), 'features', '', '', 5),
    'recoverability_atlas': ('label class', (_CLASS,), 'features', 'target:categorical', 'target', 0),
    'consensus_modules': ('gene label class', (_MAP,), 'features', 'target:categorical', 'target', 0),
    'blind_battery': ('gene label class trait', (_MAP, _REPL), 'features', '', '', 0),
    'block_ablation': ('label class', (Output('evidence_effects', 'evidence_family', 'matched held-out prediction loss when removed; not gene calls'),), 'features groups', 'target:categorical', 'target', 0),
    'feature_knn': ('gene label class', (_CALL,), 'features', 'target:categorical', 'target', 0),
    'map_neighbours': ('gene label class', (_MAP, _CALL), 'features', 'target:categorical', 'target', 0),
    'cluster_guilt': ('gene label class', (_MAP, _CALL), 'features', 'target:categorical', 'target', 0),
    'label_outliers': ('gene label class', (_RANK,), 'features', 'target:categorical', 'target', 0),
    'layer_propagation': ('gene label class', (_CALL,), 'layer', 'target:categorical', 'target', 0),
    'layer_vote': ('gene label class', (_CALL,), 'features', 'target:categorical', 'target', 0),
    'physical_partners': ('gene label class', (_CALL,), 'physical', 'target:categorical', 'target', 0),
    'structural_homology': ('gene label class', (_CALL,), 'structural layer', 'target:categorical', 'target', 0),
    'multiplex_modules': ('gene label class', (_MAP, _CALL), 'multiplex', 'target:categorical', 'target', 0),
    'link_prediction': ('gene pair', (_PAIR,), 'layer', '', '', 0),
    'attention_correction': ('gene pair label class', (_PAIR,), 'literature layer', '', '', 0),
    'unwritten_links': ('gene pair', (_PAIR,), 'measured literature', '', '', 0),
    'supervised_classifier': ('gene label class', (_CALL,), 'features groups', 'target:categorical', 'target', 0),
    'positive_unlabeled': ('gene_set', (_RANK,), 'features', '', '', 5),
    'trait_regression': ('gene trait', (_VALUE, _RANK), 'features groups', 'target:numeric', 'target', 0),
    'masked_imputation': ('gene trait', (Output('numeric_estimates', 'gene', 'imputed percentile in transformed matrix; not native assay units'),), 'features', 'column:numeric', 'column', 0),
    'condition_shift': ('gene trait', (Output('numeric_estimates', 'gene', 'rank-scale condition residual after baseline adjustment; not native assay units'),), 'features groups', 'condition:numeric baseline:numeric', 'condition', 0),
    'set_enrichment': ('gene_set', (_RANK, Output('feature_enrichment', 'feature', 'feature association and multiplicity-adjusted q; not probability of a gene label')), 'features', '', '', 3),
    'seed_expansion': ('gene_set', (_RANK,), 'features_or_graph', '', '', 2),
    'split_clusters': ('gene label class', (_MAP, _REPL), 'features', 'a:categorical b:any', 'a', 0),
    'conjunctions': ('gene label class', (_MAP, _REPL), 'features', 'a:categorical b:categorical', 'a', 0),
    'paralog_divergence': ('gene pair label class', (_PAIR,), 'features groups', '', '', 0),
    'ortholog_transfer': ('gene label class trait', (_VALUE, _CALL), 'groups other', 'target:any source:other', 'target', 0),
    'stratum_focus': ('gene label class', (_CALL,), 'features', 'target:categorical', 'target', 0),
    'triangulation': ('gene label class', (_CALL,), 'features', 'target:categorical', 'target', 0),
    'understudied_first': ('gene label class', (_CALL,), 'features attention', 'target:categorical', 'target', 0),
    'neighbour_space': ('gene gene_set pair', (_PAIR,), 'graph', '', '', 0),
    'network_training': ('gene pair', (_PAIR,), 'layer', '', '', 0),
    'conformal_calls': ('gene label class', (_CALL, Output('label_sets', 'gene', 'split-conformal prediction set; test coverage and size separately')), 'features groups', 'target:categorical', 'target', 0),
    'graph_convolution': ('gene label class', (_CALL,), 'features measured', 'target:categorical', 'target', 0),
    'random_forest': ('gene label class', (_CALL, Output('evidence_effects', 'feature', 'held-out permutation contribution; not individual component accuracy')), 'features', 'target:categorical', 'target', 0),
    'stacking': ('gene label class', (_CALL,), 'features groups', 'target:categorical', 'target', 0),
    'conformal_values': ('gene trait', (_VALUE, Output('numeric_intervals', 'gene', 'split-conformal interval in target units; coverage and width need outer testing')), 'features groups', 'target:numeric', 'target', 0),
}

_DIRECT = frozenset({'knn', 'network_vote', 'random_walk_restart', 'logistic_regression',
                     'gradient_boosting', 'ridge', 'pu_bagging', 'soft_impute',
                     'random_forest', 'stacking', 'split_conformal'})
_ABLATION = frozenset({'umap', 'hdbscan', 'co_association', 'knn_graph', 'greedy_modularity',
                       'louvain', 'tm_score', 'chance_weighting', 'agreement', 'sgc',
                       'spectral_embedding', 'attention_residual', 'support_count'})
_CONTROL = frozenset({'settings_walk', 'map_quality', 'triadic_closure', 'degree_matched',
                      'profile_correlation', 'orthogroup_mapping', 'hypergeometric',
                      'rank_sum', 'chi_square', 'kruskal_wallis', 'bh_fdr', 'grouped_cv',
                      'ablation', 'bimodal_split', 'residualisation', 'permutation_importance'})


def technique_validation(key):
    """Declare the component's test role, without borrowing pipeline accuracy."""
    if key not in _DIRECT | _ABLATION | _CONTROL:
        raise ValueError('Every technique requires an explicit reviewed validation role')
    if key in _DIRECT:
        roles = ('direct_prediction', 'pipeline_ablation')
        test = 'Held-out task truth against matched baseline; isolate component contribution on the same frozen cohort.'
    elif key in _ABLATION:
        roles = ('pipeline_ablation', 'known_truth_null_control')
        test = 'Paired pipeline with/without component on frozen biological truth; planted structure and matched null checks.'
    else:
        roles = ('known_truth_null_control',)
        test = 'Known truth and null fixtures for ' + T.TECHNIQUES[key].name + '; no biological accuracy attribution from a composite.'
    if key == 'split_conformal':
        test = 'Disjoint train/calibration/test: achieved coverage AND set/interval size; compare point-model contribution.'
    elif key in {'settings_walk', 'chance_weighting', 'grouped_cv', 'stacking'}:
        test += ' Deliberate held-out leakage must be refused; settings and ensemble fitting use inner roles only.'
    return TechniqueValidation(key, roles, test)


def catalog():
    """Return all checked declarations; undeclared new catalogue entries fail closed."""
    strategies = S.catalog()
    if {s.key for s in strategies} != set(_SPECS):
        raise ValueError('Every catalogue key requires an explicit capability declaration')
    if {t for s in strategies for t in s.techniques} != set(T.TECHNIQUES):
        raise ValueError('Every current technique requires a catalogue validation role')
    result = []
    for strategy in strategies:
        kinds, outputs, required, columns, target, seeds = _SPECS[strategy.key]
        parameters = tuple(p.name for p in strategy.params)
        column_parameters = tuple(tuple(item.split(':')) for item in columns.split())
        if set(p for p, _ in column_parameters) - set(parameters):
            raise ValueError('Capability references an absent catalogue parameter')
        tasks = (SC.T_VALUES, SC.T_LABEL) if strategy.key == 'ortholog_transfer' else (strategy.task,)
        bases = ['evidence_provenance', 'permitted_features']
        if set(required.split()) & {'layer', 'graph', 'measured', 'physical', 'structural', 'multiplex', 'literature'}:
            bases.append('permitted_graph_operators')
        if any(t in strategy.techniques for t in ('umap', 'spectral_embedding')):
            bases.append('representations_for_declared_fit_role')
        bases += ['target_and_query_specific_models_or_results', 'held_out_outputs', 'scorecards']
        result.append(Capability(strategy.key, tuple(kinds.split()), outputs, tuple(required.split()),
            column_parameters, target, seeds, parameters, tasks,
            tuple(technique_validation(t) for t in strategy.techniques), tuple(bases),
            ('organism', 'entity_universe', 'context', 'source_hashes', 'code_hashes',
             'query', 'target', 'source_exclusions', 'settings', 'fit_split', 'fit_role')))
    return tuple(result)


def get(key):
    """Get a validated capability or refuse an unknown strategy key."""
    return next((cap for cap in catalog() if cap.strategy == key), None) or _unknown(key)


def _unknown(key):
    raise KeyError('No declared capability for ' + str(key))


@dataclass(frozen=True)
class EvidenceState:
    """Explicit adapter inputs; callers must inspect data rather than infer availability.

    Columns are (name, categorical/numeric) and layers are (name, source kind):
    measured, physical, structural or literature. This is the declared evidence
    available after query-specific exclusions. Counts do not claim benchmark truth.
    """

    organism: str
    entity_ids: tuple[str, ...]
    columns: tuple[tuple[str, str], ...] = ()
    layers: tuple[tuple[str, str], ...] = ()
    permitted_feature_count: int = 0
    orthogroups: bool = False
    attention: bool = False
    other_organism: str = ''
    other_columns: tuple[str, ...] = ()
    other_orthogroups: bool = False
    class_values: tuple[tuple[str, tuple[str, ...]], ...] = ()
    contexts: tuple[BiologicalContext, ...] = (BiologicalContext(),)
    admitted_gene_space: bool = False
    entity_kind: str = 'gene'

    def __post_init__(self):
        if self.organism not in O.SPACES and self.organism not in O.HOST_TABLES:
            raise ValueError('Evidence needs an explicitly registered organism')
        if len(set(self.entity_ids)) != len(self.entity_ids) or not self.entity_ids:
            raise ValueError('Evidence needs a unique nonempty entity universe')
        if self.entity_kind not in {'gene', 'protein'} or type(self.permitted_feature_count) is not int or self.permitted_feature_count < 0:
            raise ValueError('Invalid evidence unit or permitted feature count')
        if any(type(value) is not bool for value in (self.orthogroups, self.attention, self.other_orthogroups, self.admitted_gene_space)):
            raise ValueError('Evidence availability flags must be explicit booleans')
        for values, kinds in ((self.columns, {'categorical', 'numeric'}),
                              (self.layers, {'measured', 'physical', 'structural', 'literature'})):
            if len({v[0] for v in values}) != len(values) or any(not name or kind not in kinds for name, kind in values):
                raise ValueError('Evidence names must be unique with explicit supported kinds')
        if not self.contexts or not all(isinstance(c, BiologicalContext) for c in self.contexts):
            raise ValueError('Explicit available contexts are required')
        if self.other_organism and (self.other_organism == self.organism or self.other_organism not in O.SPACES):
            raise ValueError('Orthology requires a distinct registered gene organism')
        if self.other_columns and not self.other_organism:
            raise ValueError('Other columns need their explicit organism')
        columns = dict(self.columns)
        if len(dict(self.class_values)) != len(self.class_values) or any(
                columns.get(name) != 'categorical' or not values or len(set(values)) != len(values)
                for name, values in self.class_values):
            raise ValueError('Class values need a categorical column and unique exact values')


@dataclass(frozen=True)
class AdapterPlan:
    """A deterministic answer or refusal for one query/strategy pair.

    Applicable means the declared inputs support a plan; it never means a model
    is fitted, confidence calibrated, or an independent biological test passed.
    """

    strategy: str
    query_json: str
    applicable: bool
    reason: str
    detail: str
    settings: tuple[tuple[str, object], ...] = ()
    outputs: tuple[Output, ...] = ()
    benchmark_ids: tuple[str, ...] = ()
    validation_gaps: tuple[str, ...] = ()
    schema_version: int = SCHEMA_VERSION


def _checked_settings(strategy, settings):
    params = {p.name: p for p in strategy.params}
    if set(settings) - set(params):
        raise ValueError('Unknown strategy parameter')
    for name, value in settings.items():
        p = params[name]
        if value is None and p.optional:
            continue
        if p.kind in {'int', 'float'}:
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError('Invalid numerical parameter: ' + name)
            if p.kind == 'int' and not isinstance(value, int):
                raise ValueError('Integer parameter required: ' + name)
            if p.hi > p.lo and not p.lo <= value <= p.hi:
                raise ValueError('Parameter outside declared bounds: ' + name)
        elif not isinstance(value, str):
            raise ValueError('Textual parameter required: ' + name)
        elif not value.strip() and not p.optional:
            raise ValueError('Nonempty parameter required: ' + name)
        if p.kind == 'choice' and not callable(p.choices) and value not in p.choices:
            raise ValueError('Unknown parameter choice: ' + name)


def resolve(query, key, evidence, *, settings=None, benchmarks=()):
    """Resolve one plan, preserving task meaning and explicit input/truth gaps.

    Supplied evidence must already apply the source-exclusion closure. The later
    fitting/artifact adapters must enforce split guards; this interface never
    invokes a runner, guesses a default target, or converts support to probability.
    """
    if not isinstance(query, Query) or not isinstance(evidence, EvidenceState):
        raise TypeError('Typed query and explicit evidence state are required')
    cap, strategy = get(key), S.get(key)
    supplied = dict(settings or {})
    def refused(reason, detail):
        return AdapterPlan(key, query.to_json(), False, reason, detail)
    if query.organism != evidence.organism:
        return refused('organism_mismatch', 'Query and evidence organism differ')
    if query.organism not in cap.organisms:
        return refused('organism_adapter_unavailable', 'Current catalogue adapters cover the two admitted parasite gene spaces; host gene adapters remain pending')
    if evidence.entity_kind != 'gene' or not evidence.admitted_gene_space:
        return refused('gene_space_unavailable', 'Protein references are not admitted gene spaces')
    if query.output not in {'inferences', 'scorecard'}:
        return refused('output_adapter_unavailable', 'Evidence, agreement and meta-inference use their own later adapters')
    if query.kind not in cap.query_kinds:
        return refused('query_kind_unsupported', 'Supported kinds: ' + ', '.join(cap.query_kinds))
    if any(e.organism != query.organism or e.kind != 'gene' for e in query.entities):
        return refused('endpoint_unit_or_organism_unsupported', 'This adapter requires genes in one organism; cross-organism bridges are separate')
    unknown = sorted({e.identifier for e in query.entities} - set(evidence.entity_ids))
    if unknown:
        return refused('entity_not_in_universe', ', '.join(unknown))
    if query.context not in evidence.contexts:
        return refused('context_unavailable', 'Requested context is not explicitly represented by these inputs')
    columns = dict(evidence.columns)
    if query.target and query.target not in columns:
        return refused('column_unavailable', query.target)
    if query.kind in {'label', 'class', 'trait'}:
        expected = 'numeric' if query.kind == 'trait' else 'categorical'
        if columns[query.target] != expected:
            return refused('target_kind_mismatch', 'This query requires a ' + expected + ' target')
    if query.kind == 'class':
        values = dict(evidence.class_values).get(query.target)
        if values is None:
            return refused('class_values_unavailable', 'Exact target class membership has not been inspected')
        if set(query.values) - set(values):
            return refused('class_value_unavailable', 'Requested values are absent from the target class universe')
    if cap.target_parameter and query.target:
        if cap.target_parameter in supplied and supplied[cap.target_parameter] != query.target:
            return refused('target_conflict', 'Query and supplied settings name different targets')
        supplied[cap.target_parameter] = query.target
    if query.kind == 'gene_set' or (key == 'neighbour_space' and query.kind == 'gene'):
        ids = ' '.join(e.identifier for e in query.entities)
        if 'genes' in supplied and supplied['genes'] != ids:
            return refused('seed_conflict', 'Seed parameter must preserve the exact canonical query entities')
        supplied['genes'] = ids
    if cap.minimum_seeds and len(query.entities) < cap.minimum_seeds:
        return refused('insufficient_seeds', f'At least {cap.minimum_seeds} distinct query seeds required')
    try:
        _checked_settings(strategy, supplied)
    except ValueError as error:
        return refused('invalid_settings', str(error))
    columns, layers = dict(evidence.columns), dict(evidence.layers)
    for parameter, kind in cap.column_parameters:
        column = supplied.get(parameter)
        available = evidence.other_columns if kind == 'other' else columns
        if not column:
            return refused('column_required', 'Specify ' + parameter + ' explicitly; no target is guessed')
        if column not in available:
            return refused('column_unavailable', parameter + ': ' + str(column))
        if kind not in {'any', 'other'} and columns[column] != kind:
            return refused('target_kind_mismatch', parameter + ' needs a ' + kind + ' column')
    if query.kind == 'class' and query.target not in columns:
        return refused('column_unavailable', query.target)
    flags = {
        'features': evidence.permitted_feature_count > 0,
        'features_or_graph': evidence.permitted_feature_count > 0 or bool(layers),
        'groups': evidence.orthogroups,
        'graph': bool(layers),
        'multiplex': len(layers) >= 2,
        'physical': 'physical' in layers.values(),
        'structural': 'structural' in layers.values(),
        'literature': 'literature' in layers.values(),
        'measured': bool(set(layers.values()) & {'measured', 'physical'}),
        'attention': evidence.attention,
        'other': bool(evidence.other_organism and evidence.other_orthogroups),
        'layer': supplied.get('layer') in layers,
    }
    missing = [r for r in cap.required if not flags[r]]
    if missing:
        return refused('required_evidence_unavailable', ', '.join(missing))
    if 'layer' in supplied and supplied['layer'] is not None and supplied['layer'] not in layers:
        return refused('layer_unavailable', str(supplied['layer']))
    if key == 'structural_homology' and layers.get(supplied.get('layer')) != 'structural':
        return refused('layer_kind_mismatch', 'Structural transfer requires the declared structural layer')
    if key == 'attention_correction' and layers.get(supplied.get('layer')) != 'literature':
        return refused('layer_kind_mismatch', 'Attention correction requires the declared literature layer')
    if key == 'condition_shift' and supplied['condition'] == supplied['baseline']:
        return refused('distinct_columns_required', 'A condition cannot serve as its own baseline')
    if key == 'seed_expansion' and supplied.get('mode') == 'networks' and not flags['measured']:
        return refused('required_evidence_unavailable', 'Network-only mode requires a measured network')
    outputs = cap.outputs
    if key == 'ortholog_transfer':
        outputs = (_VALUE,) if columns[supplied['target']] == 'numeric' else (_CALL,)
    applicable_tasks = (SC.T_VALUES,) if outputs == (_VALUE,) else (SC.T_LABEL,) if key == 'ortholog_transfer' else cap.benchmark_tasks
    relevant = tuple(b for b in benchmarks if b.organism == query.organism and b.task in applicable_tasks and (not query.target or b.target == query.target))
    gaps = ['Independent biological benchmark unavailable' if not any(b.status == 'admitted' for b in relevant) else 'Benchmark applicability still requires context and frozen split review',
            'Technique component tests are declared, not measured by this interface',
            'Parameters not supplied retain unresolved defaults until the context adapter binds them']
    return AdapterPlan(key, query.to_json(), True, 'applicable', 'Declared inputs support these output meanings; no inference executed',
        tuple(sorted(supplied.items())), outputs, tuple(b.benchmark_id for b in relevant), tuple(gaps))
