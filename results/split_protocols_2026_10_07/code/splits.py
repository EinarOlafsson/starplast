"""Frozen nested hold-outs and explicit fit/source guards for benchmark adapters.

Manifests describe access permissions; adapters must call the guards and record
their actual inputs. Constructing a manifest does not certify an old self-test.
Unknown homology and unavailable pair negatives are not resolved by partitioning.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from . import organisms as O

ROLES = ("train", "tune", "calibration", "test")
GROUP_KINDS = {"homology", "entity", "dataset", "context", "time"}
FIT_ACCESS = {"imputation": {"train"}, "scaling": {"train"}, "representation": {"train"},
              "feature_selection": {"train"}, "model_fit": {"train"},
              "setting_selection": {"train", "tune"}, "refit": {"train", "tune"},
              "calibration": {"calibration"}}


@dataclass(frozen=True)
class Assignment:
    """One explicit entity, protected group and assigned partition."""

    entity: str
    group: str
    role: str

    def __post_init__(self):
        if not isinstance(self.entity, str) or not self.entity or not isinstance(self.group, str) or not self.group:
            raise ValueError("Entity and protected group must be explicit nonempty strings")
        if self.role not in ROLES:
            raise ValueError("Unknown hold-out role")


@dataclass(frozen=True)
class SplitManifest:
    """A nested partition with a protected-group boundary and declared feature access."""

    organism: str
    benchmark_id: str
    group_kind: str
    assignments: tuple[Assignment, ...]
    seed: int
    feature_access: str = "inductive"
    protocol_version: int = 1

    def __post_init__(self):
        if self.organism not in O.SPACES and self.organism not in O.HOST_TABLES:
            raise ValueError("Unknown split organism")
        if not self.benchmark_id or self.group_kind not in GROUP_KINDS:
            raise ValueError("Declare benchmark identity and protected grouping")
        if type(self.seed) is not int or self.protocol_version != 1:
            raise ValueError("Unsupported split seed or version")
        if self.feature_access not in {"inductive", "transductive"}:
            raise ValueError("Declare inductive or transductive access")
        if len({row.entity for row in self.assignments}) != len(self.assignments):
            raise ValueError("Entity appears in more than one partition")
        if {row.role for row in self.assignments} != set(ROLES):
            raise ValueError("Each training/tuning/calibration/test partition must be nonempty")
        roles = {}
        for row in self.assignments:
            if row.group in roles and roles[row.group] != row.role:
                raise ValueError("Protected group crosses a hold-out boundary")
            roles[row.group] = row.role

    @property
    def identity(self) -> str:
        """Content identity including role assignment, grouping and feature-access regime."""
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    def entities(self, role) -> tuple[str, ...]:
        """Entities in a named partition, in the frozen input order."""
        if role not in ROLES:
            raise ValueError("Unknown partition")
        return tuple(row.entity for row in self.assignments if row.role == role)

    def guard_fit(self, stage, entities):
        """Reject forbidden entity access during fitting, selection or calibration."""
        if stage not in FIT_ACCESS:
            raise ValueError("Unknown fit stage; no implicit permission")
        roles = {row.entity: row.role for row in self.assignments}
        for entity in entities:
            if entity not in roles or roles[entity] not in FIT_ACCESS[stage]:
                raise ValueError("Forbidden " + stage + " access to " + str(entity))

    def guard_joint_features(self, entities):
        """Joint feature/graph construction needs explicit transductive hold-out access.

        Pointwise prediction reads held-out features without fitting; this guard
        concerns joint construction, not that legitimate inference operation.
        Learned preprocessing still uses ``guard_fit`` and training entities.
        """
        permitted = {row.entity for row in self.assignments if row.role == "train" or self.feature_access == "transductive"}
        if set(entities) - permitted:
            raise ValueError("Unlabelled hold-out feature access is not permitted by this regime")


def make_split(entities, groups, *, organism, benchmark_id, group_kind="homology", seed=0,
               fractions=(0.55, 0.15, 0.15, 0.15), feature_access="inductive", ordered=False):
    """Split entire explicit groups, using deterministic seeded or chronological order.

    Fractions apply to group counts, not to independent biological sample counts.
    For time hold-outs, callers supply sorted ISO time keys; equal times remain
    together. Dataset/context partitions protect those groups, not homology.
    """
    entities, groups = tuple(entities), tuple(groups)
    if len(entities) != len(groups) or len(set(entities)) != len(entities):
        raise ValueError("Supply one protected group per unique entity")
    for entity, group in zip(entities, groups):
        Assignment(entity, group, "train")
    fractions = np.asarray(fractions, dtype=float)
    if len(fractions) != 4 or not np.isfinite(fractions).all() or (fractions <= 0).any() or not np.isclose(fractions.sum(), 1):
        raise ValueError("Four positive role fractions must sum to one")
    unique = sorted(set(groups))
    if len(unique) < 4:
        raise ValueError("Fewer than four independent protected groups; nested split is unavailable")
    if group_kind == "time" and not ordered:
        raise ValueError("Time groups require chronological ordered=True")
    if ordered and group_kind != "time":
        raise ValueError("Chronological ordering requires explicit time groups")
    if ordered:
        from datetime import datetime
        try:
            parsed = [(datetime.fromisoformat(group), group) for group in unique]
            if len({time for time, _group in parsed}) != len(parsed):
                raise ValueError("Equivalent time keys must be normalized into one group")
            unique = [group for _time, group in sorted(parsed)]
        except (ValueError, TypeError) as error:
            raise ValueError("Time groups must be comparable ISO timestamps") from error
    else:
        unique = np.random.default_rng(seed).permutation(unique).tolist()
    # Each role gets a group before the remaining groups are apportioned.
    desired = fractions * (len(unique) - 4)
    counts = np.floor(desired).astype(int) + 1
    for index in np.argsort(-(desired - np.floor(desired)), kind="stable")[:len(unique) - int(counts.sum())]:
        counts[index] += 1
    group_roles, start = {}, 0
    for role, count in zip(ROLES, counts):
        for group in unique[start:start + int(count)]:
            group_roles[group] = role
        start += int(count)
    return SplitManifest(organism, benchmark_id, group_kind,
                         tuple(Assignment(entity, group, group_roles[group]) for entity, group in zip(entities, groups)),
                         seed, feature_access)


@dataclass(frozen=True)
class ExclusionManifest:
    """The exact target-family column and derived-layer closure for a benchmark."""

    benchmark_id: str
    targets: tuple[str, ...]
    columns: tuple[str, ...]
    layers: tuple[str, ...]
    fit_entities: tuple[str, ...]
    split_identity: str
    scope: str = "target_family"

    def __post_init__(self):
        if not self.benchmark_id or not self.targets or not set(self.targets) <= set(self.columns):
            raise ValueError("Exclusions must include the explicitly named target")
        if self.scope != "target_family":
            raise ValueError("Benchmark exclusions require target-family scope")
        if not self.fit_entities or len(set(self.fit_entities)) != len(self.fit_entities) or len(self.split_identity) != 64:
            raise ValueError("Pin the training cohort and split used for exclusion selection")

    def guard_inputs(self, *, columns=(), layers=(), benchmark_id=None):
        """Reject target, sibling-source and derived-layer contamination before fitting."""
        if benchmark_id is not None and benchmark_id != self.benchmark_id:
            raise ValueError("Exclusion manifest belongs to a different benchmark")
        leaks = (set(columns) & set(self.columns)) | (set(layers) & set(self.layers))
        if leaks:
            raise ValueError("Target/source contamination: " + ", ".join(sorted(leaks)))


def make_exclusions(nodes, targets, *, benchmark_id, split, entity_column="gene_id"):
    """Compute closure using training data only, rejecting held-out target access.

    Declared families use the full column schema; any empirical association or
    derivative detection may inspect only the supplied training rows. Final-test
    labels cannot influence feature/source exclusion selection.
    """
    from . import search
    from .strategies import DERIVED_LAYERS

    if benchmark_id != split.benchmark_id or entity_column not in nodes:
        raise ValueError("Training table and benchmark must match the frozen split")
    ids = nodes[entity_column].tolist()
    split.guard_fit("feature_selection", ids)
    if isinstance(targets, str):
        targets = (targets,)
    columns, layers = set(targets), set()
    for target in targets:
        if target not in nodes:
            raise ValueError("Named target is not in the declared table")
        columns.update(search.excluded_for(nodes, target, scope="target_family"))
        layers.update(search.excluded_layers(nodes, target, scope="target_family"))
    changed = True
    while changed:
        new = {layer for layer, parents in DERIVED_LAYERS.items() if set(parents) & layers}
        changed = bool(new - layers)
        layers.update(new)
    return ExclusionManifest(benchmark_id, tuple(targets), tuple(sorted(columns)), tuple(sorted(layers)), tuple(ids), split.identity)


@dataclass(frozen=True)
class PairPartition:
    """Positive/control pair partition with explicit endpoint regime and dropped pairs."""

    regime: str
    split: SplitManifest
    pairs: tuple[tuple[str, str, str], ...]
    excluded_cross_partition: tuple[tuple[str, str], ...] = ()

    def __post_init__(self):
        if self.regime not in {"known_node", "cold_node"}:
            raise ValueError("Declare known-node edge completion or cold-node generalization")
        seen = set()
        roles = {row.entity: row.role for row in self.split.assignments}
        train_nodes = {node for left, right, role in self.pairs if role == "train" for node in (left, right)}
        for left, right, role in self.pairs:
            key = tuple(sorted((left, right)))
            if left == right or key in seen or role not in ROLES:
                raise ValueError("Pairs must be unique undirected nonself relationships with explicit roles")
            seen.add(key)
            if self.regime == "known_node" and roles.get(pair_id(left, right)) != role:
                raise ValueError("Known-node pair role differs from its frozen edge partition")
            if self.regime == "cold_node" and (roles.get(left) != role or roles.get(right) != role):
                raise ValueError("Cold-node pair has an endpoint in another partition")
            if self.regime == "known_node" and role != "train" and not {left, right} <= train_nodes:
                raise ValueError("Known-node held-out endpoint was never observed in training")
        if self.regime == "cold_node":
            for left, right in self.excluded_cross_partition:
                if left not in roles or right not in roles or roles[left] == roles[right]:
                    raise ValueError("Dropped cold-node pair is not a cross-partition pair")


def pair_id(left, right):
    """An unambiguous undirected pair key without guessing endpoint species."""
    if not isinstance(left, str) or not left or not isinstance(right, str) or not right or left == right:
        raise ValueError("Pairs require explicit distinct nonempty endpoint identifiers")
    return json.dumps(sorted((left, right)), separators=(",", ":"))


def known_node_pairs(pairs, *, organism, benchmark_id, seed=0):
    """Reserve endpoint-covering training edges before splitting the remaining edges.

    All held-out endpoints therefore occur in training. These pairs must come
    from the caller's declared observed/control universe; missing edges are never
    generated as negatives. Sparse graphs may have no usable nested hold-out.
    """
    pairs = tuple(pairs)
    keys = [pair_id(left, right) for left, right in pairs]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate undirected pair")
    order = np.random.default_rng(seed).permutation(len(pairs))
    covered, mandatory = set(), set()
    for index in order:
        left, right = pairs[index]
        if left not in covered or right not in covered:
            mandatory.add(int(index))
            covered.update((left, right))
    remaining = [key for index, key in enumerate(keys) if index not in mandatory]
    split = make_split(remaining, remaining, organism=organism, benchmark_id=benchmark_id,
                       group_kind="entity", seed=seed, feature_access="transductive")
    roles = {row.entity: row.role for row in split.assignments}
    assignments = tuple(Assignment(key, key, "train" if index in mandatory else roles[key]) for index, key in enumerate(keys))
    full = SplitManifest(organism, benchmark_id, "entity", assignments, seed, "transductive")
    return PairPartition("known_node", full, tuple((left, right, row.role) for (left, right), row in zip(pairs, assignments)))


def cold_node_pairs(pairs, node_split):
    """Keep only same-partition edges; both test endpoints are unseen in training."""
    if node_split.group_kind not in {"homology", "entity"}:
        raise ValueError("Cold-node generalization requires node/entity or homology partitions")
    roles = {row.entity: row.role for row in node_split.assignments}
    kept, dropped, seen = [], [], set()
    for left, right in pairs:
        key = tuple(sorted((left, right)))
        if left == right or key in seen or left not in roles or right not in roles:
            raise ValueError("Pairs require unique nonself endpoints from the frozen node universe")
        seen.add(key)
        if roles[left] == roles[right]:
            kept.append((left, right, roles[left]))
        else:
            dropped.append((left, right))
    return PairPartition("cold_node", node_split, tuple(kept), tuple(dropped))


def write_split(path, split, exclusions=None):
    """Write immutable, content-identified split and matching exclusion manifests."""
    if exclusions is not None and exclusions.benchmark_id != split.benchmark_id:
        raise ValueError("Split and exclusions must name the same benchmark")
    if exclusions is not None:
        if exclusions.split_identity != split.identity:
            raise ValueError("Exclusions were selected under a different split")
        split.guard_fit("feature_selection", exclusions.fit_entities)
    payload = {"identity": split.identity, "split": asdict(split),
               "exclusions": asdict(exclusions) if exclusions else None}
    payload["content_sha256"] = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(payload, stream, indent=2)
        stream.write("\n")


def read_split(path):
    """Reconstruct partitions, checking content identity and protected-group separation."""
    payload = json.loads(Path(path).read_text())
    digest = payload.pop("content_sha256")
    if hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest() != digest:
        raise ValueError("Split/exclusion content identity changed")
    row = payload["split"]
    row["assignments"] = tuple(Assignment(**item) for item in row["assignments"])
    split = SplitManifest(**row)
    if split.identity != payload["identity"]:
        raise ValueError("Split identity changed")
    excluded = payload["exclusions"]
    exclusions = ExclusionManifest(excluded["benchmark_id"], tuple(excluded["targets"]),
                                   tuple(excluded["columns"]), tuple(excluded["layers"]), tuple(excluded["fit_entities"]),
                                   excluded["split_identity"], excluded["scope"]) if excluded else None
    if exclusions and exclusions.benchmark_id != split.benchmark_id:
        raise ValueError("Split/exclusion benchmark mismatch")
    if exclusions:
        if exclusions.split_identity != split.identity:
            raise ValueError("Exclusion selection used a different split")
        split.guard_fit("feature_selection", exclusions.fit_entities)
    return split, exclusions
