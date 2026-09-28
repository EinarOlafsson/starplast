"""About this test: four plain fields per strategy, the box under the bars on each strategy card.

Every strategy states, in a few sentences a bench biologist can read without the Guide:

    does       what it does -- the question it answers and how;
    evaluated  how it is evaluated -- what its self-test hides, what it is compared with, and the
               one number the verdict rests on;
    failure    what failure looks like, and the realistic reasons for it;
    success    what success looks like, and what new information it then gives.

Each entry was written from the strategy's own explanation, walkthrough and test description in
`strategy_catalog`, `strategy_graph` and `strategy_learning`, and claims nothing they do not. A test
requires all four fields, non-empty, for every registered strategy, so a new strategy cannot ship
without saying how it fails.
"""
from __future__ import annotations

import json
import os

#: The four fields, in the order the card shows them, with their headings.
FIELDS = (("does", "What it does"),
          ("evaluated", "How it is evaluated"),
          ("failure", "What failure looks like, and why"),
          ("success", "What success looks like, and why"))

EXPLAINERS: dict = {
    'holdout_search': {
        'does':
            ('It hides one label, such as compartment, plus anything that restates it. It then '
             'builds many maps of the remaining measurements with different settings, clusters '
             'each, and keeps the map whose clusters best isolate that label. Unlabeled genes '
             "in each label's best cluster become candidates."),
        'evaluated':
            ('A quarter of the label is hidden, whole gene families at a time. The map and '
             "each label's cluster are picked on visible genes, then scored by F1 (a 0-to-1 "
             'match score) on the hidden genes. The same clusters are rescored 100 times with '
             'hidden labels shuffled; F1 must beat that by 0.05.'),
        'failure':
            ("The hidden genes' F1 sits at the shuffled level: the chosen clusters hold the "
             'label no better than chance. Likely reasons are that the label is not written '
             'into these measurements, or that the best map was simply the luckiest of many '
             'and does not carry over to new genes.'),
        'success':
            ('Hidden genes land in the cluster chosen for their label clearly more often than '
             'shuffled labels would. The label is then a real feature of expression, fitness '
             'or modification data, and the unlabeled genes in that cluster are leads, each '
             "tagged with its cluster's F1."),
    },
    'geneset_hunt': {
        'does':
            ('You paste a gene list, such as screen hits or a complex. It walks many maps and '
             'clusterings and finds the single cluster that best captures your list, balancing '
             'how much of the cluster is on the list against how much of the list is in it. '
             "The cluster's other members are candidates."),
        'evaluated':
            ('With 20 or more genes, 30% of your list is hidden and the best cluster is chosen '
             'on the rest. The score is F1, a 0-to-1 match between the hidden genes and that '
             "cluster's other members. Twenty random lists of the same size run the same walk; "
             'yours must beat their 95th percentile by 0.05.'),
        'failure':
            ('The hidden genes score no better than random lists, which reach about 0.04. Your '
             'list is then not a unit this data can see: its genes do not behave alike in '
             'these measurements. Across many maps, something always collects part of any '
             'list, so a good-looking cluster alone proves nothing.'),
        'success':
            ('The hidden members fall back into the cluster chosen from the rest, far above '
             'random lists. Your list is then a coherent biological unit in this data. The '
             'genes that cluster with it, nearest the cluster center first, are a ranked '
             'shortlist of new members to test.'),
    },
    'recoverability_atlas': {
        'does':
            ('It builds one map with a label hidden, then asks for each category, such as each '
             "compartment, whether a labeled gene's map neighbors share its category. The "
             'result is an atlas of which categories the measurements encode and which they do '
             'not.'),
        'evaluated':
            ('30% of the label is hidden and the atlas is built from visible genes only. '
             'Hidden genes are then scored by their visible neighbors. The number is the mean '
             "AUROC (0.5 is random, 1 is perfect) for the atlas's top-half categories, versus "
             '20 shuffled-label runs, needing a 0.05 margin.'),
        'failure':
            ('The top-ranked categories score no better on hidden genes than with shuffled '
             'labels. Then the atlas does not generalize: the categories it called encoded '
             'were a fluke of the visible genes, or the label lives only in the experiment '
             'that defined it and leaves no trace elsewhere.'),
        'success':
            ('The categories the atlas ranks highest are also the ones best recovered for '
             'hidden genes. You learn which distinctions leave a real trace in expression, '
             'fitness, modification and structure, so you can take those to prediction '
             'strategies and skip searches that cannot work.'),
    },
    'consensus_modules': {
        'does':
            ('It builds and clusters many maps with different settings, then counts how often '
             'each pair of genes lands in the same cluster. Groups that stay together in at '
             'least the chosen share of maps become modules. No label is used to build them, '
             'so they reflect the data, not one lucky setting.'),
        'evaluated':
            ('Modules are built from a small walk on 1,500 genes without the label. A quarter '
             "of the label is hidden, and each label's best module is picked on visible genes. "
             'The number is the F1 match of hidden genes to that module, versus 100 label '
             'shuffles; it must beat them by 0.05.'),
        'failure':
            ('Hidden genes fit their chosen module no better than shuffled labels do. A very '
             'large module scores the same either way, so it earns nothing. Likely causes are '
             'modules too coarse to separate categories, or a label not reflected in structure '
             'that survives across maps.'),
        'success':
            ('Hidden genes sit in the module chosen for their label well above chance. That '
             'means stable modules track real biology. A stable module with no dominant known '
             'label is a candidate new unit, such as an unrecognized complex or pathway, worth '
             'following up.'),
    },
    'blind_battery': {
        'does':
            ('It builds a map from one kind of evidence only, such as transcription, tuned for '
             'clean clusters with no label in view. Then it tests every measurement the map '
             'never saw against those clusters, asking which other biology, such as fitness or '
             'compartment, the clusters separate.'),
        'evaluated':
            ('Associations are found on a random half of up to 2,000 genes and rechecked on '
             'the other half. The number is the share of findings that replicate. It is '
             "compared with 10 runs where the second half's clusters are shuffled; it must "
             'beat that by 0.2, with at least three findings.'),
        'failure':
            ('Few findings replicate, no more than with shuffled clusters, or fewer than three '
             'are found. The significant results were then large-sample artifacts: with many '
             'genes, tiny differences look significant. Or the chosen evidence simply does not '
             'couple to the other measurements.'),
        'success':
            ('Most associations found on one half hold on the other, evidence that two kinds '
             'of biology are coupled. For example, co-expressed clusters that differ in '
             'fitness suggest transcription programs track essentiality. Another family may '
             'encode something different.'),
    },
    'block_ablation': {
        'does':
            ('For a label, it asks which kind of evidence carries it. Each family of '
             'measurements, or each single dataset, predicts the label alone by a '
             '15-nearest-neighbor vote. It also measures how much the full combination loses '
             'when that evidence is removed.'),
        'evaluated':
            ('A quarter of the label is hidden. Each kind of evidence is ranked by accuracy on '
             'visible labels, and the top one predicts the hidden labels. That hidden accuracy '
             'is compared with the hidden accuracy of every kind of evidence, as if picked at '
             'random; it must top their 80th percentile.'),
        'failure':
            ('The top-ranked evidence does no better on hidden genes than a randomly chosen '
             'kind. The ranking is then not real: several kinds of evidence may carry the '
             'label about equally, or none of them carries it well, so the order on visible '
             'genes was noise.'),
        'success':
            ('The evidence ranked first also predicts hidden labels best. You learn which '
             'experiments encode this biology, which are redundant (strong alone, no loss when '
             'removed) and which are unique (costly to remove). That guides a focused map and '
             'where the next experiment should go.'),
    },
    'feature_knn': {
        'does':
            ('For each unlabeled gene, it finds the most similar labeled genes across all '
             'permitted measurements at once. Those neighbors vote, closer ones counting more. '
             'The gene is called with the winning label if it holds enough of the vote. No map '
             'or clustering sits in between.'),
        'evaluated':
            ('A quarter of the label is hidden, whole gene families at a time. Each hidden '
             'gene is called by its nearest visible genes. The number is the share of hidden '
             'genes called correctly, with abstentions as misses. It must beat 10 '
             'shuffled-label runs by 0.05.'),
        'failure':
            ('The share called correctly sits at the shuffled level. With hundreds of '
             'measurements, many missing and filled with the median, the nearest genes may '
             'share gaps in the data rather than biology. Or the label is not encoded in these '
             'measurements at all.'),
        'success':
            ('Hidden genes are called correctly well above chance. This sets the baseline '
             'other strategies must beat. Each call on an unlabeled gene is transparent: you '
             'can trace it to the neighbors that made it, and a stricter vote threshold gives '
             'fewer, safer calls.'),
    },
    'map_neighbours': {
        'does':
            ('It builds one map with the label withheld. Each unlabeled gene is then called by '
             'the distance-weighted vote of its nearest labeled neighbors on the map. It is '
             'what you do by eye when you see a grey point inside a colored cloud, made '
             'systematic.'),
        'evaluated':
            ('On up to 2,500 genes the map places, a quarter of the label is hidden, whole '
             'gene families at a time. Hidden genes are called by their nearest visible map '
             'neighbors. The number is the share called correctly, compared with 20 '
             'shuffled-label runs, needing a 0.05 margin.'),
        'failure':
            ('Correct calls sit at the shuffled level. The map may have torn the category '
             'apart when squeezing hundreds of measurements into three dimensions, or the '
             'label is not encoded in the data. Genes left out of a sampled map are never '
             'called.'),
        'success':
            ('Hidden genes are called correctly well above chance. Comparing with strategy 07 '
             'on the same label tells you whether the map cleans up noise or loses signal. The '
             'calls label unlabeled genes, and you can check each in context on the colored '
             'map.'),
    },
    'cluster_guilt': {
        'does':
            ('It clusters a map built without the label, then tests every cluster against '
             'every label for enrichment, asking whether a label appears more often than '
             'chance. Unlabeled members of clusters that are both significant and enriched by '
             'the chosen lift are called with that label.'),
        'evaluated':
            ('A quarter of the label is hidden on up to 2,500 mapped genes, and enrichment is '
             'recomputed from visible labels. The number is precision: how often a call on a '
             'hidden gene is right, since the method abstains by design. It must beat 20 '
             'shuffled-label runs by 0.1.'),
        'failure':
            ('Calls on hidden genes are right no more often than with shuffled labels. An '
             'enriched cluster can still have most labeled members elsewhere, so its calls can '
             'be mostly wrong. Clusters too small to hold several labeled genes also give weak '
             'enrichment.'),
        'success':
            ('Calls on hidden genes are right clearly more often than chance. A cluster can '
             'then be named by its enrichment, such as an apicoplast-flavored cluster, even if '
             'it is not pure. Its unlabeled members become candidates, each backed by a '
             'q-value and a lift.'),
    },
    'label_outliers': {
        'does':
            ('It turns the question around and asks which existing labels look wrong. For each '
             'labeled gene it checks how many of its nearest measurement neighbors and network '
             'partners share its label. Genes whose company mostly carries another label rank '
             'as most surprising.'),
        'evaluated':
            ('It deliberately swaps 5% of labels, at least ten, to a wrong class, then '
             'computes surprise. The number is the AUROC (0.5 is random, 1 is perfect) for the '
             'swapped genes rising to the top. It is compared with 20 random gene sets of the '
             'same size and must beat them by 0.1.'),
        'failure':
            ('Swapped genes rank no higher than random genes. Neighbors then do not agree '
             'enough on the label to expose a wrong one, because the label is weakly encoded '
             'in the measurements or too few genes have labeled network partners.'),
        'success':
            ('Deliberately swapped labels rise to the top of the list. On real data, the top '
             'genes are likely annotation errors, proteins with two locations or functions, or '
             'genes with unusual measurements. Each comes with the label its neighbors carry '
             'instead.'),
    },
    'layer_propagation': {
        'does':
            ('It seeds each label on the genes that carry it and lets it spread along the '
             'edges of one measured network, such as crosslinks or co-fitness, by a random '
             'walk. Hubs are down-weighted. Each reached gene is called by the label that '
             'arrives most strongly.'),
        'evaluated':
            ('A quarter of the label is hidden, whole gene families at a time, and only '
             'visible genes seed the spread. The number is the share of hidden genes called '
             'correctly, with unreached genes as misses. It is compared with 10 shuffled-label '
             'runs and must beat them by 0.05.'),
        'failure':
            ('Correct calls sit at the shuffled level. The chosen network does not link genes '
             'that share this label, or it reaches too few genes, since unreached genes count '
             'as misses. A layer built from the label itself is refused, as it would be '
             'circular.'),
        'success':
            ('Hidden genes are called correctly well above chance, so this kind of '
             'relationship carries the label. Testing each layer builds a table of which '
             "networks encode which biology. The best layer's calls label the unlabeled genes "
             'it reaches.'),
    },
    'layer_vote': {
        'does':
            ("Every permitted network and the measurement neighbors vote on each gene's label. "
             "Each source's vote is weighted by how far it beat its own chance level on known "
             'labels. The verdict reaches more genes than any single network and leans on '
             'sources that have earned trust.'),
        'evaluated':
            ('A quarter of the label is hidden. Weights are learned on an inner holdout of '
             'visible labels only, then hidden genes are called. The number is the share '
             'called correctly, compared with 5 shuffled-label runs where weights are '
             'relearned; it must beat them by 0.05.'),
        'failure':
            ('Correct calls sit at the shuffled level, and source weights are near zero. No '
             'source knew more than chance about this label, so the label is invisible in '
             'these networks and measurements, however many edges they have.'),
        'success':
            ('Hidden genes are called correctly well above chance. The source weights table '
             'ranks each kind of evidence by how much it knows about this label, net of '
             "chance. The calls label unlabeled genes, with support as the winner's share of "
             'the weighted vote.'),
    },
    'physical_partners': {
        'does':
            ('It calls each protein by the labels of its measured physical partners, from '
             'crosslinks and replicated pulldowns, weighted by how often each contact was '
             'seen. Crosslinked proteins were close together in the parasite, so they share a '
             'compartment. Every call lists its partners.'),
        'evaluated':
            ('Only genes with at least one physical partner are scored. A quarter of the label '
             'is hidden, whole gene families at a time, and hidden genes are called by their '
             "visible partners' weighted vote. The share called correctly must beat 20 "
             'shuffled-label runs by 0.05.'),
        'failure':
            ('Correct calls sit at the shuffled level. Few hidden genes may have labeled '
             'partners, since the interactomes cover a few thousand proteins at most. The '
             "crosslinker's chemistry also favors lysines and soluble proteins, so some "
             'proteins are poorly sampled.'),
        'success':
            ('Hidden proteins are placed correctly well above chance. Within its reach, this '
             'is expected to be the most precise localization evidence available. Each call '
             'names the partners, labels and contact counts behind it, so you can check it by '
             'hand.'),
    },
    'structural_homology': {
        'does':
            ("It asks what a protein's 3D fold says about its enzyme class, even when its "
             'sequence matches nothing annotated. Predicted structures are compared, and each '
             'unannotated protein takes a vote of its annotated look-alikes, weighted by how '
             'similar the folds are (TM-score).'),
        'evaluated':
            ('Among proteins with a structural neighbor, a quarter of the annotations are '
             'hidden, whole gene families at a time, and called back from their neighbors. The '
             'same test is run 20 times on shuffled annotations. The verdict rests on the '
             'share of hidden proteins called correctly.'),
        'failure':
            ('The correct-call rate sits at the shuffled level. Likely reasons: the chosen '
             'depth is too fine, since a fold rarely knows the exact substrate at EC levels 3 '
             'or 4, or too few proteins have structural neighbors, or the layer chosen is not '
             'the structural one.'),
        'success':
            ('The correct-call rate clears the shuffled runs, which says fold carries this '
             'annotation at this depth. You get calls for hypothetical proteins, listed first, '
             "with the neighbors' models to check. Testing a deeper level shows how far that "
             'trust extends.'),
    },
    'multiplex_modules': {
        'does':
            ('It finds groups of genes that form communities in several kinds of relationship '
             'at once, such as co-expression, co-fitness and crosslinks. Communities are found '
             'in each network separately, and two genes share a module only when enough of '
             'those networks put them together.'),
        'evaluated':
            ('A quarter of a chosen label is hidden. On the visible genes, the module that '
             'best matches each label is picked, and the test asks how well it holds that '
             "label's hidden genes (an F1 score). The null reshuffles the hidden genes' labels "
             '100 times.'),
        'failure':
            ('The F1 score is no better than with shuffled labels. The modules may reflect a '
             'bias the networks share, such as favoring abundant proteins, rather than '
             'biology. Or the chosen label simply does not follow these communities, or the '
             'resolution makes modules too big or small.'),
        'success':
            ('The hidden genes land in the module chosen for their label, so the agreement '
             'between networks carries real signal. The modules become candidate complexes and '
             'pathways. A large module with no dominant label is a set of genes worth reading '
             'one by one.'),
    },
    'link_prediction': {
        'does':
            ('It predicts physical contacts that a crosslinking or pulldown experiment missed. '
             'A model learns what separates measured contacts from other pairs: shared '
             'partners, support in other networks, and similar measurements. It then scores '
             'every unobserved pair with any support.'),
        'evaluated':
            ("A fifth of the layer's contacts are hidden and the model is retrained without "
             'them. It must rank the hidden contacts above fresh pairs of genes with a similar '
             'number of partners. The null is 5 models fed scrambled gene evidence; the number '
             'is the AUROC (ranking accuracy).'),
        'failure':
            ('The hidden contacts rank no better than with scrambled evidence. That means the '
             "other layers do not know about this interactome's contacts beyond each gene's "
             'popularity, or the layer has too few edges to learn from. Derived and literature '
             'layers are refused as targets.'),
        'success':
            ('The hidden contacts outrank matched non-contacts, so the other evidence can '
             'point to real missed interactions. You get a ranked list of likely contacts with '
             'the layers backing each. A link from a labeled to an unlabeled protein also '
             'hints at its location.'),
    },
    'attention_correction': {
        'does':
            ('It ranks gene pairs by how often papers mention them together beyond what each '
             "gene's fame predicts. Each co-mention count is replaced by its excess over what "
             "the two genes' own publication totals would give, so famous genes stop "
             'dominating the top.'),
        'evaluated':
            ('The label is never used to build the literature layer. Among pairs with both '
             'genes labeled, the top corrected pairs are taken, and the number is the share '
             'sharing a label. It is compared with 20 random sets of co-mentioned pairs; the '
             'raw ranking is shown beside it.'),
        'failure':
            ('The top corrected pairs share a label no more often than random co-mentioned '
             "pairs. Then the literature's links do not track this label, or too few "
             'co-mentioned pairs have both genes labeled for the top to be a real top.'),
        'success':
            ('The corrected top pairs share a compartment or complex more often than chance, '
             'so the literature links them for a reason. Pairs with a high score but no shared '
             "label are the literature's hypotheses that the current labels do not yet "
             'explain.'),
    },
    'unwritten_links': {
        'does':
            ('It counts, for every gene pair, how many independent measurement networks link '
             'it, then lists the well-supported pairs that no paper has ever mentioned '
             'together. Paralogs, shared domains and the compartment layer are left out, since '
             'those links are already written down.'),
        'evaluated':
            ('The literature serves only as the answer key. The test asks whether pairs with '
             'more measurement support are co-mentioned more often than random pairs, scored '
             'as an AUROC (ranking accuracy). The null repeats this 10 times with gene '
             'identities scrambled.'),
        'failure':
            ('The AUROC sits at the scrambled level: measurement support does not predict what '
             'biologists write about on this table. Then the unwritten pairs are not a '
             'forecast of future papers, and the list is just agreement between experiments, '
             'not a knowledge gap.'),
        'success':
            ('Better-supported pairs are co-mentioned more often, so the unwritten ones at the '
             'top are the most likely to be written next. Each is a concrete hypothesis with '
             'its supporting layers named. A pair of one understudied and one well-studied '
             'gene is cheapest to follow up.'),
    },
    'supervised_classifier': {
        'does':
            ('It trains a model on the labeled genes to learn which measurements recognize '
             'each class, then assigns a label with a probability to every unlabeled gene. It '
             'is balanced so the most common class does not win by default, and it lists the '
             'measurements each class relies on.'),
        'evaluated':
            ('A quarter of the label is hidden, whole gene families at a time, so a paralog '
             'cannot give the answer away. The model calls the hidden genes. The number is the '
             'share called correctly, compared with the chance rate for a guesser with the '
             'same mix of answers.'),
        'failure':
            ('The correct-call rate is near chance. The label may not be encoded in the '
             'permitted measurements, or the penalty setting makes the model too simple or '
             'lets it fit noise. A class recognized mainly by missing data is another warning '
             'sign.'),
        'success':
            ('Hidden genes are called well above chance, so the measurements do carry the '
             'label. You get probability-ranked calls for unlabeled genes and a readable list '
             'of what defines each class. Raising the probability threshold leaves fewer but '
             'surer calls.'),
    },
    'positive_unlabeled': {
        'does':
            ('It starts from a list of genes that ARE something, with no list of genes that '
             'are not. Many small models each compare your list with a small random draw of '
             'other genes, and every other gene is scored only by models that did not train on '
             'it, then ranked.'),
        'evaluated':
            ('With 20 or more genes, 30% of your list is hidden and the rest trains the '
             'models; smaller lists are tested on a known category of similar size. The number '
             'is the AUROC: how well hidden members outrank every other gene. The null is 10 '
             'random lists of the same size.'),
        'failure':
            ('Hidden members rank no better than for random lists. The measurements may not '
             'capture what your list shares, or the list may be too mixed to have a common '
             'signature. If the list came from a column, failing to name it would instead '
             'inflate the result.'),
        'success':
            ('Hidden members rise to the top, so the models learned what sets your list apart, '
             'even from weak signals spread over many measurements. The candidates are genes '
             'that look like list members; take the top ones to strategy 24 to see what they '
             'share.'),
    },
    'trait_regression': {
        'does':
            ('It predicts one numeric measurement, such as a fitness score, from all other '
             'permitted measurements. By default, other measurements of the same kind are left '
             'out. Each gene is predicted by a model that never saw it or its paralog.'),
        'evaluated':
            ('A fifth of the measured values are hidden and the model is trained on the rest. '
             'The number is the rank correlation between predicted and hidden values. The null '
             'is 3 models trained on the same values shuffled among genes.'),
        'failure':
            ('The correlation is no better than with shuffled values. The measurement then has '
             'no signature in the other data, so predictions for unmeasured genes mean little '
             'and the surprising genes may be noise. Including screens of the same kind would '
             'mainly show agreement.'),
        'success':
            ('The hidden values are predicted, so the trait has a signature in the other data. '
             'Unmeasured genes get an estimated value, trusted as far as the correlation says. '
             'Genes far from their prediction are candidates for unusual function, or '
             'artifacts worth rechecking.'),
    },
    'masked_imputation': {
        'does':
            ('It estimates missing measurements from the rest of the table using a model of a '
             'few shared patterns (a low-rank model). It first checks, column by column, how '
             'well each measurement can be rebuilt, and fills only where that works; elsewhere '
             'a gap stays a gap.'),
        'evaluated':
            ("A tenth of every column's measured values is hidden and the table is completed. "
             'Each column gets a reliability, the rank correlation on its hidden values; the '
             'verdict uses the median over columns. The null is 3 tables with each column '
             'shuffled on its own.'),
        'failure':
            ('The median reliability is at the shuffled level. The measurements are then not '
             'correlated enough to predict each other, or the rank is off: too few patterns '
             'blur programs together, too many fit noise. Missing values should stay unknown '
             'rather than be filled.'),
        'success':
            ('Columns are rebuilt far better than shuffled ones, so the measurements share '
             'structure. Columns with high reliability can be filled for genes never measured, '
             'which beats filling with the median. The reliability list also says which '
             'columns must stay gaps.'),
    },
    'condition_shift': {
        'does':
            ('It finds genes that matter more or less in one condition, such as the mouse, '
             'than a baseline like fibroblasts predicts. The condition screen is corrected for '
             'the baseline, and the leftover shift is then predicted from other measurements, '
             'with both screens withheld.'),
        'evaluated':
            ('A fifth of the genes measured in both screens are hidden, and a model trained on '
             'the rest predicts their shift. The number is the rank correlation between '
             'predicted and actual shifts. The null is 3 models trained on shuffled shifts.'),
        'failure':
            ('The correlation stays at the shuffled level. Then the condition-specific shift '
             'is idiosyncratic or noise, with no signature in the other data. The ranked list '
             'of shifted genes remains, but predictions for genes never screened in the '
             'condition should not be used.'),
        'success':
            ('Hidden shifts are predicted, so the condition-specific need has a signature, '
             'such as secretion or host-facing location. Strongly negative genes are '
             'candidates for host interaction or in vivo nutrition, and unscreened genes can '
             'be ranked by that signature.'),
    },
    'set_enrichment': {
        'does':
            ('It describes what your gene list has in common, testing every category, '
             'measurement and network in the table at once with one multiple-testing '
             'correction. The significant features form a profile, which is then used to rank '
             'every other gene.'),
        'evaluated':
            ('The profile is built from 60% of the list, and the other 40% is hidden. The '
             'number is the AUROC: how well hidden members outrank every other gene on the '
             'profile. The null is 10 random lists of the same size.'),
        'failure':
            ('Hidden members score no better than random lists. The enriched features are then '
             'a description of this particular list rather than something its members truly '
             'share. Naming the column the list came from matters, or the profile just reads '
             'that column back.'),
        'success':
            ('The profile recognizes members it never saw, so it captures what the list truly '
             'shares. You get a readable profile, sorted by q-value, and a ranked set of genes '
             'that resemble the list, as candidates for membership.'),
    },
    'seed_expansion': {
        'does':
            ('It grows your gene list along the measured networks. A random walk starts at '
             'your genes, steps along edges and jumps back often; genes it visits most are '
             'close by many short paths. Layers are averaged, and hubs are down-weighted so '
             'they do not attract every walk.'),
        'evaluated':
            ('30% of your list is hidden and the walk starts from the other 70%. The number is '
             'the AUROC: how well hidden members outrank every other gene by visit score. The '
             'null is 10 random starting lists of the same size.'),
        'failure':
            ('Hidden members are visited no more than for random lists. Your genes may have '
             'few network edges, or be linked only through single long paths. The measurement '
             'graph option can reach genes no network covers.'),
        'success':
            ('The walk returns to the hidden members, so your list is tightly knit in the '
             "data. The candidates are the list's closest relatives, such as other members of "
             'a complex, pathway or secretory route, each shown with the layers linking it '
             'directly to your genes.'),
    },
    'split_clusters': {
        'does':
            ('It finds clusters that agree on one label, such as a compartment, but split '
             'cleanly on another label or measurement, such as a stage or a fitness level. '
             'Both labels are withheld from the map, so each split is a claim that a category '
             'holds two kinds of gene.'),
        'evaluated':
            ('Splits are found on a random half of the genes and checked on the other half. '
             'The number is the share of findings that reappear. The null is 20 runs with the '
             'splitting label shuffled in the second half; passing also needs at least three '
             'findings.'),
        'failure':
            ('Few findings replicate, no more than with the shuffled label, or there are fewer '
             'than three. The category may have no real internal structure, or clusters are '
             'too coarse or too small to hold a whole category with both sides. A fine-grained '
             'shared label rarely clusters.'),
        'success':
            ('Splits reappear in genes the discovery never saw, so the category truly contains '
             'two kinds of gene, and the measurement carrying the split is named. The minority '
             'side is listed, a lead on, for example, which proteins in a compartment differ '
             'by phase or essentiality.'),
    },
    'conjunctions': {
        'does':
            ('It looks for groups of genes defined by two labels at once, such as secreted AND '
             'fitness-conferring. Both labels are hidden, genes are clustered on a map of the '
             "measurements, and each cluster's joint enrichment is compared with the stronger "
             'of the two single enrichments. A ratio above 1.25 flags a combination.'),
        'evaluated':
            ('Combinations are found on a random half of the genes and checked on the other '
             'half: among genes with the first label, is the second label still concentrated '
             'in that cluster? The number that decides it is the share of findings that '
             'replicate, compared with 20 runs where the second label is shuffled.'),
        'failure':
            ('The replicating share sits near the shuffled level, or fewer than three '
             'combinations are found. Often the cluster was simply enriched for one label, so '
             'the second adds nothing within it. Clusters may also be too coarse to isolate a '
             'combination, or the two labels are not encoded in these measurements.'),
        'success':
            ('Several combinations replicate well above the shuffled level on genes the '
             'discovery never saw. That means a real kind of gene exists that neither label '
             "shows alone. The cluster's unlabeled members become candidates for the "
             'combination.'),
    },
    'paralog_divergence': {
        'does':
            ('It asks which duplicated genes (paralogs) behave differently across the '
             'measurements, hinting that one copy took a new job. Divergence is one minus the '
             "correlation of the two genes' profiles over the measurements both have. Pairs "
             'are ranked from most to least diverged.'),
        'evaluated':
            ('Pairs where both genes carry a compartment label are used, with that label '
             'hidden from the profiles. The test asks whether pairs in different compartments '
             'are more diverged than pairs in the same one. The deciding number is that '
             'ranking score (AUROC), versus 20 random reassignments of divergence.'),
        'failure':
            ('The score sits at the random-reassignment level: divergence does not track a '
             'change of compartment. Pairs may share too few measurements, making divergence '
             'noisy. Or copies in different compartments may still be measured alike, so the '
             'profiles cannot tell them apart.'),
        'success':
            ('Pairs in different compartments are clearly more diverged than pairs in the same '
             'one. Divergence then carries meaning, so a highly diverged pair with no '
             'localization is evidence that one copy changed its place or job. It gives a '
             'ranked shortlist of paralogs to compare side by side.'),
    },
    'ortholog_transfer': {
        'does':
            ('It carries a measurement or label from the other parasite onto this one through '
             'shared orthogroups (gene families across species). The source is summarized per '
             'orthogroup, its relation to the target is learned on genes measured in both, and '
             'that is applied to genes measured only on the other side.'),
        'evaluated':
            ('A quarter of the genes measured in both species are hidden and the relation is '
             'learned on the rest. For a measurement, the deciding number is the rank '
             'correlation between transferred and hidden values; for a label, correct calls. '
             'Both are compared with runs where ortholog values go to random genes.'),
        'failure':
            ('The transfer does no better than orthologs dealt to random genes. The two '
             'species may use the gene differently, so the measurement is not conserved. Too '
             'few genes may be measured in both to learn the relation. Lineage-specific genes '
             'have no ortholog at all and are never reached.'),
        'success':
            ('Transferred values track the hidden ones well above the shuffled level. The '
             'ortholog then acts as a second, independent measurement of the gene. A '
             'Plasmodium knockout screen, for example, becomes a prediction for untested '
             "Toxoplasma genes, listed in 'transferred'."),
    },
    'stratum_focus': {
        'does':
            ('It calls labels for one hard group of genes, such as lineage-specific, '
             'hypothetical or understudied ones, from their nearest neighbors in the '
             'measurements. All conservation and orthology columns are removed first, so a '
             'call cannot just reflect how conserved a gene is.'),
        'evaluated':
            ('A quarter of the label is hidden by whole orthogroups, and only hidden genes '
             'inside the chosen group are scored. The deciding number is correct calls per '
             'hidden gene of that group. It is compared with 10 runs on shuffled labels, and '
             'shown beside the accuracy for the rest of the proteome.'),
        'failure':
            ('Accuracy within the group sits at the shuffled level, even if the whole-proteome '
             'number looked good. These genes are where methods are weakest: few labeled '
             'neighbors resemble them, and the label may not be encoded in their measurements '
             'once conservation is removed.'),
        'success':
            ('Calls in the group beat shuffled labels by a clear margin. The error rate then '
             'belongs to the genes being predicted, not borrowed from well-studied ones. '
             "Unlabeled genes in the group get calls that can be trusted at that group's own "
             'rate.'),
    },
    'triangulation': {
        'does':
            ('It runs three methods built on different evidence: measurement neighbors, a '
             'trained classifier and network partners. A gene is called only where enough of '
             'them give the same label, two of three by default. Their failures are largely '
             'independent, so agreement should be right more often.'),
        'evaluated':
            ('A quarter of the label is hidden; each method trains on the visible genes. The '
             'deciding number is precision, the share of agreed calls on hidden genes that are '
             'correct. It is compared with agreement between the same methods trained on '
             "shuffled labels, and each method's own precision is shown."),
        'failure':
            ('Agreed calls are no more precise than chance agreement on shuffled labels. The '
             'methods may share a failure, so they agree on wrong answers. Or agreement '
             'reaches very few genes, since a call needs every method to speak, and few genes '
             'may have network partners.'),
        'success':
            ('Agreed calls are clearly more precise than chance and than single methods. The '
             'result is a smaller, surer set of labels for unlabeled genes. Raising the '
             'agreement to three of three gives the smallest and most certain set.'),
    },
    'understudied_first': {
        'does':
            ('It makes the agreed calls of strategy 31 only for genes with no focal or '
             'substantive paper. Candidates are ranked by how many methods agree times how '
             'little the gene has been written about. The top of the list is where the data '
             'can say something new.'),
        'evaluated':
            ('A quarter of the label is hidden, and agreed calls are scored only on hidden '
             'understudied genes. The deciding number is precision, the share of those calls '
             'that are right. It is compared with 5 runs on shuffled labels, never borrowing '
             'trust from well-studied genes.'),
        'failure':
            ('Precision on understudied genes sits at the shuffled level. The measurements '
             'were often designed around well-studied genes, so understudied ones may be '
             'poorly covered. The methods may also rarely agree on them, leaving too few calls '
             'to trust.'),
        'success':
            ('Agreed calls on understudied genes are clearly more precise than chance. Each '
             'top candidate is then a first hypothesis about a gene with no literature. The '
             'precision shown beside it is the only evidence for that call there is.'),
    },
    'neighbour_space': {
        'does':
            ('It merges every permitted network layer and the measurement table into one graph '
             "and lists a gene's nearest neighbors. Each pair gets one probability, with the "
             'sources that support it. Nothing is written back into a layer, and each pair is '
             'marked as measured or inferred.'),
        'evaluated':
            ("One layer's edges are hidden by whole orthogroups and that layer is dropped from "
             'the inputs. The deciding number is how well hidden edges rank above non-pairs of '
             'similar connectivity (AUROC). It is compared with 10 rewirings that keep each '
             "gene's edge count but scramble which genes are linked."),
        'failure':
            ('The score does not beat the rewired networks. The graph then mostly knows which '
             'genes are well connected, not who pairs with whom. The other layers and '
             'measurements may carry little about the hidden layer, especially if it is '
             'sparse.'),
        'success':
            ('Hidden edges are found clearly above the rewired level. Then the combined graph '
             "knows real pairings, not just busy genes. A gene's neighbor list, with the "
             'evidence behind each entry, can suggest partners that no single network records.'),
    },
    'network_training': {
        'does':
            ('It trains one model per network layer, each blind to its own layer, and ranks '
             'gene pairs that the evidence implies but no layer records. Each gap comes with a '
             'calibrated probability and the sources behind it. An interpretable baseline '
             'ranks unless a learned embedding clearly beats it.'),
        'evaluated':
            ("The chosen layer's edges are hidden by orthogroup and removed from the inputs. "
             'The deciding number is how well hidden edges rank above non-pairs of similar '
             'connectivity (AUROC), against 10 rewirings that keep edge counts. The easier '
             'random-pair score and their gap are shown beside it.'),
        'failure':
            ('A score near 0.5, or at the rewired level, means the model learned only which '
             'genes are well connected. Crosslinks show a subtler failure: complexes are '
             'recovered, but rewired pairs inside a complex score almost as well, so single '
             'missing pairs cannot be named.'),
        'success':
            ("Hidden edges clear the rewired bar. The 'gaps' table then gives a ranked "
             'shortlist of pairs no experiment recorded, each with its supporting evidence. '
             'These are claims to test, for example with strategy 13 or 16, not measured '
             'edges.'),
    },
    'conformal_calls': {
        'does':
            ("It turns a classifier's scores into a set of possible labels for each gene, "
             'promised to hold the true label for at least 1 - alpha of new genes (90% by '
             'default). A one-label set is a call with that promise behind it. A two-label set '
             'shows the data cannot tell those compartments apart.'),
        'evaluated':
            ('A quarter of the label is hidden by orthogroup, and the model trains and '
             'calibrates on the rest. The deciding number is set efficiency: how far the sets '
             'narrow the possible labels. It is compared with 10 shuffled-label runs, and '
             'coverage on hidden genes is checked against the promise.'),
        'failure':
            ('Sets are as wide as on shuffled labels, so the measurements barely narrow the '
             'answer. Coverage can also fall short if new genes do not resemble the '
             'calibration genes. A rare compartment may be under-covered unless per-class '
             'thresholds are used.'),
        'success':
            ('Sets are clearly narrower than chance while coverage meets the promise. Many '
             'unlabeled genes then get single-label calls with a stated error rate. The '
             'multi-label sets are findings too, pointing to compartments these data cannot '
             'separate.'),
    },
    'graph_convolution': {
        'does':
            ("It averages each gene's measurements over its network partners, one and two "
             "steps out, and trains a logistic regression on the gene's own profile beside its "
             "neighborhood's. Only measurements travel along edges, never labels, so hidden "
             'labels cannot leak through neighbors.'),
        'evaluated':
            ('A quarter of the label is hidden by whole orthogroups and the model trains on '
             'the visible genes. The deciding number is the share of hidden genes called '
             'correctly. It is compared with the chance level for the same mix of predicted '
             'and true classes.'),
        'failure':
            ('Accuracy sits at the chance level for those class mixes. The label may not be '
             'encoded in the measurements or the networks. Genes with no edges fall back to '
             'their own profile, so sparse networks add little, while two steps on a dense '
             'network blur most genes together.'),
        'success':
            ("Calls clearly beat chance, and ideally beat strategy 19, which uses the gene's "
             "own profile only. The 'where the model looks' table then shows which "
             "compartments the networks encode better than the gene's own measurements. "
             'Unlabeled genes get calls.'),
    },
    'random_forest': {
        'does':
            ('It trains a random forest, hundreds of voting decision trees, on every permitted '
             'measurement to call a label. Trees can capture thresholds and combinations, such '
             'as high expression AND a signal peptide. Measurements are ranked by how much '
             'accuracy drops when each is shuffled.'),
        'evaluated':
            ('A quarter of the label is hidden by whole orthogroups and the forest trains on '
             'the rest. The deciding number is the share of hidden genes called correctly. It '
             'is compared with the chance level for the same mix of predicted and true '
             'classes.'),
        'failure':
            ('Accuracy sits at chance: the measurements do not define the label. If the forest '
             "merely ties strategy 19, the signal is linear and the simpler model's weights "
             'explain it better. Noisy labels can also let trees memorize training genes when '
             'the leaf size is small.'),
        'success':
            ('Calls clearly beat chance and ideally beat strategy 19. Unlabeled genes then get '
             "calls, and 'what defines the label' names the measurements that truly matter. "
             'Shuffling any of them costs accuracy on genes the forest never saw.'),
    },
    'stacking': {
        'does':
            ('It learns how much to trust each kind of evidence for this label: measurement '
             'neighbors, a logistic regression and network partners. Each predicts labeled '
             'genes out of fold, and a second model learns which to believe for which class. '
             'It then calls unlabeled genes.'),
        'evaluated':
            ('A quarter of the label is hidden by whole orthogroups; the base predictions use '
             'only visible genes, so no hidden label reaches either level. The deciding number '
             'is the share of hidden genes called correctly, compared with the chance level '
             'for the same class mixes.'),
        'failure':
            ('Accuracy sits at chance, meaning none of the three evidence types carries the '
             'label. If it only ties strategies 07, 19 or 31, combining adds little for this '
             'label. Few folds give the second model noisier training data.'),
        'success':
            ("Calls clearly beat chance, and usually beat the single methods. 'Trust by "
             "evidence' then answers which experiments matter for this label, for example "
             'networks for complexes but not secreted proteins. Unlabeled genes get the '
             'combined calls.'),
    },
    'conformal_values': {
        'does':
            ('It predicts a measurement, such as a fitness score, from the other permitted '
             "measurements. Each prediction gets an interval in the measurement's own units, "
             'promised to contain the true value for at least 1 - alpha of new genes. Measured '
             'genes outside their interval are listed as surprises.'),
        'evaluated':
            ('A fifth of the measured values are hidden, and the model trains and calibrates '
             'on the rest. The deciding number is the rank correlation between predicted and '
             'hidden values, compared with the chance level. The share of hidden values inside '
             'their intervals is checked against the promise.'),
        'failure':
            ('Predictions track hidden values no better than chance, and intervals are as wide '
             "as the measurement's spread. The other measurements then know little about this "
             'one. Leaving out the same kind of screen can remove the only related evidence.'),
        'success':
            ('Predictions track hidden values and intervals are narrower than the spread. '
             'Unmeasured genes get a value worth acting on, with a stated range. Measured '
             'genes far outside their interval are surprises the rest of the data cannot '
             'explain.'),
    },
}


def explainer(key: str) -> dict:
    """The four fields for one strategy, or empty strings if it has none (the test forbids that)."""
    e = EXPLAINERS.get(key, {})
    return {f: e.get(f, "") for f, _title in FIELDS}


def first_sentence(text: str) -> str:
    """The first sentence of a field: the one line the collapsed box shows."""
    text = " ".join(str(text).split())
    for i, ch in enumerate(text):
        if ch in ".!?" and i > 20 and (i + 1 == len(text) or text[i + 1] == " "):
            return text[:i + 1]
    return text


def markdown(key: str) -> str:
    """One strategy's explainer as a Markdown list, for the docs."""
    e = explainer(key)
    return "\n".join(f"- **{title}.** {e[f]}" for f, title in FIELDS)


# --------------------------------------------------------------------------- worked examples
#: One real failure and one real success per strategy and organism, chosen from the calibration
#: sweep by `scripts/build_strategy_examples.py` (see its notebook for every choice).
EXAMPLES = os.path.join(os.path.dirname(__file__), "data", "strategy_examples.json")
_EXAMPLES: dict = {}


def examples(key: str, organism: str, path: str = EXAMPLES) -> dict:
    """{"failure": {...}, "success": {...}, "task": ...} for one strategy, or {} if none shipped."""
    if path not in _EXAMPLES:
        try:
            with open(path) as fh:
                _EXAMPLES[path] = json.load(fh)
        except (OSError, ValueError):
            _EXAMPLES[path] = {}
    return dict((_EXAMPLES[path].get(organism) or {}).get(key) or {})


def example_verdict(example: dict) -> dict:
    """An example's verdict numbers in the form `scorecard.headline_bars` takes."""
    return {"skill": example.get("skill"), "observed": example.get("observed"),
            "chance": example.get("chance")}
