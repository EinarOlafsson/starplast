"""Start here: the question tree that turns what a user HAS into the strategies worth running.

The Strategies tab is a catalogue. A catalogue answers "what is there", and someone who arrives with
a gene list from a screen and no idea which of thirty-nine named methods applies to it is not helped
by being shown all thirty-nine. This module is the other way round: it asks what the person has, what
they want to know, and ends at a short ranked list of strategies with their settings already filled
in from the answers.

It holds no Qt, so the tree and the recommendations are testable headless and can be driven from a
notebook::

    from starplast import guided, strategies
    ctx = strategies.shipped("Tg")                   # any Context
    answers = {"have": guided.HAVE_LABEL, "space": ctx.organism, "label": "compartment",
               "goal": guided.PREDICT}
    for r in guided.recommend(answers, ctx):
        print(r.number, r.title, "--", r.why)

**The tree.** One question per screen, and a screen's options may depend on every answer before it:

    have      what the person has: one gene, a gene list, their own measurement, a label, nothing
    space     which organism's table to work in (the registry's available spaces)
    gene /    the subject itself: the gene, the gene set, the label column or the numeric column,
    genes /   the last two listed with their coverage so an empty column is never chosen blind
    label /
    measure
    label     asked a second time only when a single gene is to be given a label
    goal      what they want out of it, from the goals that make sense for what they have
    baseline  a second measurement, asked only when the goal is to compare two conditions

:func:`route` is the whole branching, :func:`next_step` the one question to ask now, and an answer
dictionary with no `next_step` is finished. Every first answer reaches an end state, and every option
leads somewhere: `tests/test_guided.py` walks the tree exhaustively to check exactly that.

**The recommendations.** :func:`recommend` ranks strategies by, in order of weight:

1. whether the strategy is one of the handful this goal and this kind of subject are FOR (`KEYS`),
   best first;
2. whether its scorecard task is one the goal asks for (`GOALS[...].tasks`);
3. what the calibration sweep measured for it ON THIS SPACE (`GRADE_WEIGHT`) -- a strategy that is
   reliable for Toxoplasma and weak for Plasmodium is recommended differently in the two;
4. how long it takes, as a tie-breaker only.

A strategy whose settings cannot be filled from this table at all (no numeric column, no edge layer,
no second species) is dropped before ranking rather than recommended and then failing, and a strategy
that cannot use the subject the person chose is dropped too -- a label-calling method is not offered
for a numeric screen. :func:`caveat` says plainly when the best thing available is only weak, and
:func:`views` offers the map gallery or the star map where those answer the question better than any
strategy does.

Nothing here computes anything about the data beyond reading the columns' coverage: it decides what
to offer, and `strategies` does the work.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import scorecard as SC

# --------------------------------------------------------------------------- what the user has
HAVE_GENE = "one gene"
HAVE_SET = "a gene list"
HAVE_NUMBER = "my own measurement"
HAVE_LABEL = "a label"
HAVE_NOTHING = "nothing yet"

#: The subject step each first answer leads to -- the thing the rest of the questions are about.
SUBJECT = {HAVE_GENE: "gene", HAVE_SET: "genes", HAVE_NUMBER: "measure", HAVE_LABEL: "label",
           HAVE_NOTHING: ""}

# --------------------------------------------------------------------------- goals
MORE_LIKE = "more like mine"
PREDICT = "predict it"
EXPLAIN = "explain it"
PARTNERS = "partners"
COMPARE = "compare"
LEARNABLE = "learnable"
#: The goal of someone who chose "nothing yet": show what this table can do. Never asked, always
#: implied, so that first answer reaches an end state with one question fewer.
TOUR = "show me"


@dataclass(frozen=True)
class Option:
    """One answer to one question: what is stored, what is shown, and the tooltip that explains it."""
    value: str
    label: str
    hint: str = ""
    #: A second line under the label -- a column's coverage, a set's size. Never load-bearing.
    detail: str = ""


@dataclass(frozen=True)
class Step:
    """One screen: its question in a few plain words, and what kind of answer it takes.

    :ivar kind: ``choice`` (pick one of :func:`options`), ``gene`` (one identifier), ``genes`` (a
        list, however it is supplied) or ``column`` (a choice whose options carry coverage).
    """
    key: str
    question: str
    kind: str
    tip: str
    note: str = ""


@dataclass(frozen=True)
class Goal:
    """One thing a person can want, the scorecard tasks that answer it, and how to say why."""
    key: str
    label: str
    hint: str
    tasks: tuple
    #: The clause a recommendation's reason is built around: "<method> <reason>".
    reason: str
    #: The goal named as a thing being done: "the first choice for <noun>".
    noun: str = ""


GOALS = {
    MORE_LIKE: Goal(MORE_LIKE, "Find more genes like mine",
                    "Rank every other gene by how much it resembles the ones you gave, and say how "
                    "often that ranking has been right on gene sets whose answer was known.",
                    (SC.T_RANK, SC.T_SET),
                    "ranks every other gene by how much it looks like yours",
                    "finding more genes like yours"),
    PREDICT: Goal(PREDICT, "Predict it for the genes nobody has measured",
                  "Learn the chosen label or measurement from the genes that have one, and give a "
                  "value to the genes that have none.",
                  (SC.T_LABEL, SC.T_VALUES),
                  "gives the unmeasured genes a value and says how often that value is right",
                  "predicting it for the unmeasured genes"),
    EXPLAIN: Goal(EXPLAIN, "Explain what defines it",
                  "Say which measurements carry the thing you chose, which are redundant, and which "
                  "genes sit where they should not.",
                  (SC.T_LABEL, SC.T_RANK, SC.T_VALUES, SC.T_REPL),
                  "says which measurements carry it and which add nothing",
                  "explaining what defines it"),
    PARTNERS: Goal(PARTNERS, "Find partners or complex members",
                   "Rank the genes most likely to be physically or functionally tied to yours, with "
                   "the evidence each link rests on.",
                   (SC.T_RANK,),
                   "ranks the genes most likely to be tied to yours, with the evidence for each",
                   "finding partners"),
    COMPARE: Goal(COMPARE, "Compare two conditions",
                  "Ask which genes matter more in one condition than the other measurements say "
                  "they should -- and whether the rest of the table explains which.",
                  (SC.T_VALUES, SC.T_REPL),
                  "finds the genes that shift between your two measurements more than the rest of "
                  "the table explains",
                  "comparing two conditions"),
    LEARNABLE: Goal(LEARNABLE, "Check whether it is learnable at all",
                    "Before trusting any prediction: is this thing in the measurements at all, and "
                    "for which of its categories?",
                    (SC.T_CLUSTER, SC.T_LABEL, SC.T_REPL, SC.T_RANK),
                    "measures whether the data holds it at all, before anything is predicted from "
                    "it",
                    "checking whether it is learnable at all"),
    TOUR: Goal(TOUR, "Show me what this data can do",
               "The strategies that need nothing from you: what the measurements organise, and "
               "what they link that nobody has written down.",
               (SC.T_REPL, SC.T_RANK, SC.T_CLUSTER),
               "needs nothing from you and says what the measurements already organise",
               "a first look at what this data holds"),
}

#: The goals offered for each first answer, in the order they are shown.
GOALS_FOR = {
    HAVE_GENE: (PARTNERS, MORE_LIKE, PREDICT),
    HAVE_SET: (MORE_LIKE, PARTNERS, EXPLAIN, LEARNABLE),
    HAVE_LABEL: (PREDICT, EXPLAIN, LEARNABLE),
    HAVE_NUMBER: (PREDICT, COMPARE, EXPLAIN, LEARNABLE),
    HAVE_NOTHING: (),
}

#: The strategies each (goal, subject) pair is FOR, best first. Everything else that fits the goal's
#: task can still be recommended, below these, so a thin pair is never a dead end.
KEYS = {
    (MORE_LIKE, "gene"): ("seed_expansion", "neighbour_space", "link_prediction", "set_enrichment"),
    (MORE_LIKE, "genes"): ("set_enrichment", "seed_expansion", "positive_unlabeled",
                           "neighbour_space", "geneset_hunt", "network_training"),
    (PARTNERS, "gene"): ("link_prediction", "unwritten_links", "neighbour_space", "seed_expansion",
                         "structural_homology"),
    (PARTNERS, "genes"): ("link_prediction", "unwritten_links", "neighbour_space", "seed_expansion",
                          "multiplex_modules"),
    (PREDICT, "gene"): ("feature_knn", "layer_propagation", "triangulation",
                        "supervised_classifier", "physical_partners", "conformal_calls"),
    (PREDICT, "label"): ("supervised_classifier", "conformal_calls", "random_forest", "stacking",
                         "graph_convolution", "feature_knn", "cluster_guilt", "layer_propagation",
                         "triangulation"),
    (PREDICT, "measure"): ("trait_regression", "conformal_values", "masked_imputation"),
    (EXPLAIN, "genes"): ("set_enrichment", "geneset_hunt", "blind_battery"),
    (EXPLAIN, "label"): ("block_ablation", "recoverability_atlas", "supervised_classifier",
                         "label_outliers", "blind_battery"),
    (EXPLAIN, "measure"): ("trait_regression", "masked_imputation", "blind_battery"),
    (LEARNABLE, "genes"): ("geneset_hunt", "set_enrichment"),
    (LEARNABLE, "label"): ("recoverability_atlas", "holdout_search", "block_ablation",
                           "supervised_classifier", "blind_battery"),
    (LEARNABLE, "measure"): ("trait_regression", "masked_imputation", "conformal_values"),
    (COMPARE, "measure"): ("condition_shift", "trait_regression", "split_clusters"),
    (TOUR, ""): ("blind_battery", "recoverability_atlas", "unwritten_links", "link_prediction",
                 "consensus_modules", "label_outliers"),
}

#: What a calibration grade is worth in the ranking. A grade is measured per space, so the same
#: strategy is recommended differently for two organisms -- which is the point of measuring it.
GRADE_WEIGHT = {"reliable": 1.5, "works when tuned": 0.8, "weak": -0.4, "no skill": -3.0,
                "untestable": -1.0, "": 0.0}
#: What it is worth that a strategy takes the very genes the person chose, rather than ranking the
#: whole table and leaving them to find their own genes in it.
SUBJECT_WEIGHT = 0.7
#: Cost breaks a tie and nothing more: a faster strategy first when two are otherwise equal.
COST_WEIGHT = {"seconds": 0.15, "a minute": 0.05, "minutes": 0.0}
#: The grades a person can rely on. Below this the recommendation is still made -- there may be
#: nothing better -- but `caveat` says so in plain words.
TRUSTED = ("reliable", "works when tuned")
#: How many strategies a recommendation list holds. Short enough to read, long enough to choose.
MAX_RECOMMENDED = 6
#: The smallest list worth showing before strategies outside `KEYS` are added to fill it out.
MIN_RECOMMENDED = 3


# --------------------------------------------------------------------------- the steps
STEPS = {
    "have": Step("have", "What do you have?", "choice",
                 "Everything else follows from this answer, and you can come back and change it "
                 "from the trail at the top."),
    "space": Step("space", "Which organism?", "choice",
                  "The table the questions are asked of. Each species is measured, and its "
                  "strategies calibrated, on its own."),
    "gene": Step("gene", "Which gene?", "gene",
                 "Type an accession, a symbol or a word from the product description, and pick the "
                 "gene from the list.",
                 "search by id, symbol or product"),
    "genes": Step("genes", "Which genes?", "genes",
                  "The set the strategies start from: paste it, load a file, take the genes gated "
                  "on the map, or try an example set of a known category.",
                  "paste, load a file, use the gated genes, or try an example"),
    "label": Step("label", "Which label?", "column",
                  "The annotation the strategies work on. Each is shown with how many genes carry "
                  "it and how many classes it has, so an almost-empty column is never picked by "
                  "accident.",
                  "how many genes carry each one is shown"),
    "measure": Step("measure", "Which measurement?", "column",
                    "The number per gene the strategies work on -- a screen, a fitness score, an "
                    "abundance. Each is shown with how many genes have a value.",
                    "how many genes have a value is shown"),
    "goal": Step("goal", "What do you want to know?", "choice",
                 "The goal decides which strategies are recommended and how they are ranked."),
    "baseline": Step("baseline", "Compared with what?", "column",
                     "The second measurement: the baseline the first one is compared against, such "
                     "as growth in the dish behind growth in the animal.",
                     "the measurement the first is compared against"),
}


def answered(answers, key: str) -> bool:
    """Whether `key` holds a real answer -- an empty list or an empty string is not one."""
    value = (answers or {}).get(key)
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    return bool(len(value)) if hasattr(value, "__len__") else True


def subject(answers) -> str:
    """Which subject step this branch uses: `gene`, `genes`, `label`, `measure` or none."""
    return SUBJECT.get((answers or {}).get("have"), "")


def goal_of(answers) -> str:
    """The goal in force: the chosen one, or the tour that "nothing yet" implies."""
    if (answers or {}).get("have") == HAVE_NOTHING:
        return TOUR
    return (answers or {}).get("goal") or ""


def route(answers) -> list:
    """The steps this branch asks, in order, as far as the answers so far determine them.

    The list grows as questions are answered: an unanswered step is the last entry, because what
    comes after it is not yet known. When every step in the route is answered the route is complete
    and :func:`next_step` returns nothing.
    """
    answers = answers or {}
    chain = ["have"]
    if not answered(answers, "have"):
        return chain
    chain.append("space")
    if not answered(answers, "space"):
        return chain
    sub = subject(answers)
    if sub:
        chain.append(sub)
        if not answered(answers, sub):
            return chain
    if answers.get("have") == HAVE_NOTHING:
        return chain                                   # the tour asks nothing more
    chain.append("goal")
    if not answered(answers, "goal"):
        return chain
    # A single gene given a label needs to know WHICH label; a gene set does not, because the
    # strategies that expand a set rank genes rather than label them.
    if sub == "gene" and answers["goal"] == PREDICT:
        chain.append("label")
        if not answered(answers, "label"):
            return chain
    if answers["goal"] == COMPARE:
        chain.append("baseline")
    return chain


def next_step(answers) -> str:
    """The one question to ask now, or "" when the answers are complete and the list is due."""
    for key in route(answers):
        if not answered(answers or {}, key):
            return key
    return ""


def done(answers) -> bool:
    """Whether this branch has reached its end state: the recommendations."""
    return next_step(answers) == ""


def trail(answers) -> list:
    """(step, value, shown) for every answer so far, in the order they were given.

    The progress trail at the top of the panel, and the way back to any earlier question.
    """
    out = []
    for key in route(answers):
        if not answered(answers, key):
            break
        value = answers[key]
        out.append((key, value, shown(key, value)))
    return out


def shown(step: str, value) -> str:
    """One answer in a few words, as the trail shows it."""
    if step in ("genes",) and not isinstance(value, str):
        return f"{len(value):,} genes"
    if step == "space":
        from . import organisms
        try:
            return organisms.get(str(value)).species
        except KeyError:
            return str(value)
    if step == "goal":
        goal = GOALS.get(str(value))
        return goal.label if goal else str(value)
    text = str(value)
    if step in ("label", "measure", "baseline"):
        return text.replace("_", " ")
    if len(text) > 42:
        return text[:39] + "..."
    return text


# --------------------------------------------------------------------------- what each step offers
def options(step: str, answers, ctx=None) -> list:
    """The options of one step, given the answers before it and the table they are asked of.

    A `gene` or `genes` step has no fixed options -- the answer is an identifier or a list -- and
    returns an empty tuple, which is how the view knows to draw a search box instead of buttons.
    """
    answers = answers or {}
    if step == "have":
        return list(_HAVE_OPTIONS)
    if step == "space":
        return _space_options()
    if step == "goal":
        return [Option(g, GOALS[g].label, GOALS[g].hint)
                for g in GOALS_FOR.get(answers.get("have"), ())]
    if step == "label":
        return categorical_options(ctx)
    if step == "measure":
        return numeric_options(ctx)
    if step == "baseline":
        return [o for o in numeric_options(ctx) if o.value != answers.get("measure")]
    return []


_HAVE_OPTIONS = (
    Option(HAVE_GENE, "A single gene",
           "One accession or symbol. The questions then ask what it is tied to, what it resembles, "
           "or what label it should carry."),
    Option(HAVE_SET, "A list of genes",
           "A screen's hits, a complex, a pathway -- anything you can paste, load from a file, or "
           "draw around on the map."),
    Option(HAVE_NUMBER, "My own measurement or screen",
           "A number per gene: a fitness score, an abundance, a phenotype strength. Pick the "
           "column it is in, or import yours first under Analysis.", "a number per gene"),
    Option(HAVE_LABEL, "A label I want explained or extended",
           "An annotation some genes carry and most do not -- a localization, a phenotype class, a "
           "family.", "an annotation some genes carry"),
    Option(HAVE_NOTHING, "Nothing specific -- show me what this can do",
           "The strategies that need nothing from you: what the measurements organise on their "
           "own, and what they link that no paper mentions."),
)


def _space_options() -> list:
    """Every organism whose table is installed, as the registry names it."""
    from . import organisms
    out = []
    for code in organisms.codes(available=True):
        space = organisms.get(code)
        out.append(Option(code, space.species,
                          f"{space.species} ({space.reference}), curated by {space.database}. Its "
                          f"strategies are calibrated on its own table.", space.kind))
    return out


def categorical_options(ctx) -> list:
    """Every label that can be held out, each with its coverage: the space's own held-out labels
    first, then the rest by how many genes carry them.

    Coverage is on the option, not behind it, because the commonest way to waste an afternoon here
    is to pick a column that eleven genes carry; and the space's declared targets lead because a
    column with a value for every gene is usually a derived one, not a measured one.
    """
    if ctx is None:
        return []
    from . import organisms
    try:
        preferred = tuple(organisms.get(ctx.organism).targets)
    except KeyError:
        preferred = ()
    out = []
    for column in ctx.categorical_columns():
        truth = ctx.truth(column).dropna()
        n, classes = len(truth), truth.nunique()
        out.append(Option(column, column.replace("_", " "),
                          f"{n:,} of {ctx.n:,} genes carry this label, in {classes:,} classes. "
                          f"The commonest are {', '.join(map(str, truth.value_counts().index[:3]))}.",
                          f"{n:,} genes · {classes:,} classes"))
    order = {c: i for i, c in enumerate(preferred)}
    return sorted(out, key=lambda o: (order.get(o.value, len(order)), -_count(o.detail)))


def numeric_options(ctx) -> list:
    """Every numeric column worth predicting, the space's own screens first, then by coverage."""
    if ctx is None:
        return []
    from . import organisms
    try:
        preferred = tuple(organisms.get(ctx.organism).numbers)
    except KeyError:
        preferred = ()
    out = []
    for column in ctx.numeric_targets():
        values = ctx.values(column).dropna()
        out.append(Option(column, column.replace("_", " "),
                          f"{len(values):,} of {ctx.n:,} genes have a value, from "
                          f"{values.min():.3g} to {values.max():.3g}.",
                          f"{len(values):,} genes"))
    order = {c: i for i, c in enumerate(preferred)}
    return sorted(out, key=lambda o: (order.get(o.value, len(order)), -_count(o.detail)))


def _count(detail: str) -> int:
    """The gene count an option's detail line begins with, for sorting; 0 when it has none."""
    head = detail.split(" ")[0].replace(",", "")
    return int(head) if head.isdigit() else 0


def example_genes(ctx, answers=None) -> tuple:
    """(category, gene ids): a real gene set to try the gene-list branch on, from the default label.

    The same set `strategies.example_set` gives the Strategies tab's "example" button, so the two
    offer the same thing.
    """
    from . import strategies as S
    target = (answers or {}).get("label") or S.default_category(ctx)
    category, positions = S.example_set(ctx, target)
    return category, [str(g) for g in ctx.gene_ids[positions]]


def read_gene_file(path: str) -> list:
    """Every identifier in a gene-list file: the first column of a CSV or TSV, or a plain list.

    Qt-free, so a list can be loaded in a script exactly as the panel's "Load a file" button loads
    it, and so the reading is testable without a file dialog.
    """
    with open(path, errors="replace") as fh:
        lines = [line.strip() for line in fh if line.strip()]
    tokens = [line.replace("\t", ",").split(",")[0].strip().strip('"') for line in lines]
    return [t for t in tokens if t]


def gene_ids(answers, ctx) -> list:
    """The genes the answers name, resolved against this table; empty when none were given."""
    raw = (answers or {}).get("genes") or (answers or {}).get("gene") or []
    if not len(raw):
        return []
    positions, _missing = ctx.resolve_genes(raw)
    return [str(g) for g in ctx.gene_ids[positions]]


# --------------------------------------------------------------------------- the recommendations
@dataclass(frozen=True)
class Recommendation:
    """One recommended strategy: why it is here, and the settings the answers already decided."""
    key: str
    number: int
    title: str
    method: str
    task: str
    grade: str
    score: float
    why: str
    settings: dict = field(default_factory=dict)
    #: True when nothing better than a weak grade was available for this goal on this space.
    weak: bool = False


def usable(strategy, ctx) -> bool:
    """Whether this table can fill everything the strategy needs before a user chooses anything.

    Asked of the strategy's own defaults rather than of a list kept here: a required column, layer
    or second species that resolves to nothing on this table means the strategy has nothing to run
    on, whatever the user answers. Gene lists are exempt -- the user supplies those.
    """
    try:
        defaults = strategy.defaults(ctx)
    except Exception:
        return False
    for p in strategy.params:
        if p.optional or p.kind == "genes":
            continue
        if p.kind in ("category", "number", "column", "layer") and not defaults.get(p.name):
            return False
    return True


def fits_subject(strategy, answers, ctx) -> bool:
    """Whether the strategy can actually take the subject the person chose.

    A label-calling method offered for someone's numeric screen would silently run on a default
    label instead of on their measurement, and report a number about something they never asked
    about. That is the failure this rule exists to prevent.
    """
    sub = subject(answers)
    kinds = {p.kind for p in strategy.params}
    if sub == "measure":
        return bool(kinds & {"number"}) or strategy.task == SC.T_VALUES
    if sub in ("gene", "genes"):
        if goal_of(answers) in (MORE_LIKE, PARTNERS):
            return "genes" in kinds or strategy.task in (SC.T_RANK, SC.T_SET)
        return True
    if sub == "label":
        return bool(kinds & {"category", "column"})
    return True


def grade_of(key: str, organism: str) -> str:
    """What the calibration sweep graded this strategy on this space; "" if it never measured it."""
    from . import calibration as CAL
    return str((CAL.entry(key, organism) or {}).get("grade") or "")


def _grade_phrase(key: str, organism: str, grade: str) -> str:
    """The grade in a clause: what it means and, where measured, the skill behind it."""
    from . import calibration as CAL
    from . import organisms
    try:
        species = organisms.get(organism).species
    except KeyError:
        species = organism
    if not grade:
        return f"not calibrated on {species} yet, so test it here before believing it"
    entry = CAL.entry(key, organism) or {}
    skill = ((entry.get("default") or {}).get("skill"))
    runs = entry.get("runs") or 0
    measured = (f" (skill {skill:.2f} over {runs:,} held-out tests)"
                if isinstance(skill, (int, float)) and runs else "")
    if grade == "reliable":
        return f"reliable on {species}{measured}"
    if grade == "works when tuned":
        return f"weak at its defaults on {species} but reliable at the tuned setting{measured}"
    if grade == "weak":
        return (f"only weak on {species}{measured}: some skill on average, not enough to rely on")
    if grade == "no skill":
        return f"no better than shuffled data on {species}{measured}"
    return f"untestable on {species}: too few of its tests could be scored here"


def prefill(strategy, answers, ctx) -> dict:
    """The settings the answers decide, as `Strategy.settings` takes them.

    Only parameters the strategy actually has are filled, and only from an answer that means the
    same thing: the gene list into its gene-list parameter, the chosen label into its held-out
    label, the chosen measurement into the number it predicts, the baseline into the baseline.
    Everything else stays at the strategy's own default.
    """
    names = {p.name: p for p in strategy.params}
    out = {}
    genes = gene_ids(answers, ctx)
    if genes:
        for p in strategy.params:
            if p.kind == "genes":
                out[p.name] = list(genes)
    label = (answers or {}).get("label")
    if label and label in ctx.nodes:
        for name in ("target", "a"):
            p = names.get(name)
            if p is not None and p.kind in ("category", "column"):
                out[name] = label
    measure = (answers or {}).get("measure")
    if measure and measure in ctx.nodes:
        for name in ("target", "condition", "column"):
            p = names.get(name)
            if p is not None and p.kind in ("number", "column"):
                out[name] = measure
    baseline = (answers or {}).get("baseline")
    if baseline and baseline in ctx.nodes and "baseline" in names:
        out["baseline"] = baseline
    return out


def _score(strategy, goal: Goal, rank, grade: str, uses_subject: bool = False) -> float:
    """A strategy's place in the ranking: what it is for, what task it does, how it graded, cost.

    `uses_subject` is the one thing the goal and the grade cannot see: whether the strategy actually
    takes the genes the person chose. One that does is preferred over one that ranks the whole table
    and leaves them to look their own genes up in it.
    """
    base = (3.0 - 0.25 * rank) if rank is not None else 0.6
    if strategy.task in goal.tasks:
        base += 0.8
    if uses_subject:
        base += SUBJECT_WEIGHT
    return base + GRADE_WEIGHT.get(grade, 0.0) + COST_WEIGHT.get(strategy.cost, 0.0)


def _why(strategy, goal: Goal, answers, grade_phrase: str, ctx, rank) -> str:
    """One line: what it asks, why it is on this list, what of yours it uses, and how it graded.

    What it asks is the strategy's OWN question, whole, rather than a sentence kept here that would
    drift away from the catalogue as strategies change -- and rather than its first clause, which
    for half of them is a subordinate one and reads as an unfinished sentence.
    """
    clause = str(strategy.question).strip()
    if rank == 0:
        fit = f"The first choice for {goal.noun or goal.label.lower()}"
    elif rank is not None:
        fit = f"One of the methods for {goal.noun or goal.label.lower()}"
    else:
        fit = f"Offered because its {strategy.task} task fits {goal.noun or goal.label.lower()}"
    sub = subject(answers)
    if sub in ("gene", "genes") and any(p.kind == "genes" for p in strategy.params):
        genes = gene_ids(answers, ctx)
        fit += f", starting from your {len(genes):,} gene{'s' if len(genes) != 1 else ''}"
    elif sub == "label" and answers.get("label"):
        fit += f", on {answers['label'].replace('_', ' ')}"
    elif sub == "measure" and answers.get("measure"):
        fit += f", on {answers['measure'].replace('_', ' ')}"
    stop = "" if clause.endswith(("?", ".", "!")) else "."
    # The track record, where it covers this strategy and the label in hand: a measured rate the
    # reader can weigh, rather than only a grade.
    try:
        from . import track_record
        organism = (answers or {}).get("space") or ctx.organism
        label = answers.get("label") if sub == "label" else None
        record = track_record.record_phrase(strategy.key, organism, label)
    except Exception:
        record = ""
    record = f"; {record}" if record else ""
    return f"{clause}{stop} {fit}; {grade_phrase}{record}."


def recommend(answers, ctx, limit: int = MAX_RECOMMENDED) -> list:
    """The strategies worth running for these answers, best first, each with its reason.

    Ranked by :func:`_score`; anything this table cannot fill, or that cannot take the subject the
    person chose, is dropped first. When the handful a goal is FOR leaves fewer than
    `MIN_RECOMMENDED`, strategies whose task matches the goal fill the list rather than leaving it
    short -- so no answer is ever a dead end.
    """
    from . import strategies as S
    goal = GOALS.get(goal_of(answers))
    if goal is None:
        return []
    organism = (answers or {}).get("space") or ctx.organism
    chosen_genes = bool(gene_ids(answers, ctx))
    preferred = KEYS.get((goal.key, subject(answers)), ())
    order = {key: i for i, key in enumerate(preferred)}
    rows = []
    for strategy in S.catalog():
        rank = order.get(strategy.key)
        if rank is None and strategy.task not in goal.tasks:
            continue
        if not usable(strategy, ctx) or not fits_subject(strategy, answers, ctx):
            continue
        grade = grade_of(strategy.key, organism)
        takes = chosen_genes and any(p.kind == "genes" for p in strategy.params)
        rows.append((rank, strategy, grade, _score(strategy, goal, rank, grade, takes)))
    chosen = [r for r in rows if r[0] is not None]
    if len(chosen) < MIN_RECOMMENDED:
        chosen += [r for r in rows if r[0] is None]
    # A strategy whose own track record on this label is no better than always naming the commonest
    # class goes to the back, whatever its grade: the grade came from a different test, and the
    # record is the direct answer to "would it have known?".
    try:
        from . import track_record
        label = (answers or {}).get("label") if subject(answers) == "label" else None
        below = {r[1].key: track_record.beats_baseline(r[1].key, organism, label) is False
                 for r in chosen}
    except Exception:
        below = {}
    chosen.sort(key=lambda r: (below.get(r[1].key, False), -r[3], r[1].number))
    top = chosen[:max(1, int(limit))]
    best = max((GRADE_WEIGHT.get(g, 0.0) for _rank, _s, g, _sc in top), default=0.0)
    weak = best < GRADE_WEIGHT["works when tuned"]
    out = []
    for rank, strategy, grade, score in top:
        out.append(Recommendation(
            key=strategy.key, number=strategy.number, title=strategy.name, method=strategy.method,
            task=strategy.task, grade=grade, score=round(float(score), 3),
            why=_why(strategy, goal, answers, _grade_phrase(strategy.key, organism, grade), ctx,
                     rank),
            settings=prefill(strategy, answers, ctx), weak=weak))
    return out


def caveat(recommendations) -> str:
    """One plain sentence when the best thing available is not one to rely on, else "".

    Said rather than hidden: a person who ran the top recommendation and read a number off it
    should know beforehand that nothing better than "weak" exists here for what they asked.
    """
    if not recommendations:
        return ("Nothing here fits that combination on this table. Go back a step and change the "
                "goal, or pick a column more genes carry.")
    if not any(r.grade in TRUSTED for r in recommendations):
        best = recommendations[0]
        grade = best.grade or "not calibrated"
        return (f"None of these is reliable for this on this table: the best, {best.title}, is "
                f"{grade}. Run it if you like, but treat what it returns as a hypothesis to test, "
                f"not a result -- press Test to measure it on your own choice of label.")
    return ""


#: The other two places an answer can be better served than by any strategy: the pregenerated maps
#: and their per-label scores, and one gene's links with their provenance.
VIEW_MAPS = "maps"
VIEW_STAR = "star map"


def views(answers, ctx=None) -> list:
    """Panels, not strategies, that answer this goal better -- each with the reason it is offered."""
    goal, sub = goal_of(answers), subject(answers)
    out = []
    if goal in (LEARNABLE, EXPLAIN, TOUR) and sub in ("label", ""):
        out.append({"view": VIEW_MAPS, "title": "Map gallery",
                    "why": "Every shipped map scored against every label: which map separates the "
                           "one you chose, and how well, without running anything."})
    if sub in ("gene", "genes") or goal in (PARTNERS, TOUR):
        out.append({"view": VIEW_STAR, "title": "Star map",
                    "why": "Your gene in the middle and everything linked to it around it, each "
                           "link coloured by the measurement or strategy that claims it."})
    return out
