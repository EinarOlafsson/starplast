# Saving an annotation, with the numbers that justify it

Depends on [20_validation_completion.md](20_validation_completion.md). **Must not ship before it.**

## The workflow

A cluster is homogeneous for a category; its unlabelled members are candidates. Select the cluster,
list its members split by label / no-label, mark the unlabelled ones, persist.

`validate.candidates()` already returns exactly those rows, with `cluster_frac_category` and
`cluster_frac_contradicting` attached. What is missing is the store and the interface.

## The store

A separate annotations file. Every row carries: the gene, the proposed label, the configuration and
cluster it came from, the cluster's composition, the validation precision and recall **for that
category**, the date, and free-text reasoning.

**Never written back into the node table.** Never coloured like a measured call. This project's whole
discipline is that measurement, inference and absence never read as one another; an annotation is a
fourth thing and needs its own colour and its own file.

## Why the ordering constraint is not negotiable

Annotation is easy to build and easy to build wrongly. Shipping it before validation gives the
project a mechanism for manufacturing unvalidated claims, which is the one thing it exists to
prevent. Measured on the coarse level-of-detail clustering, annotating from those components would be
right about 1-6% of the time -- and the candidate list looks identical either way.

## Done when

A candidate can be saved and reloaded, always displays its validation numbers, and is drawn in a
colour used for nothing else.
