# Published processed review inputs for legacy derivations

NIH NLM NCBI PubMed Central Article Datasets supplied exact article media URLs
and MD5 checksums. The executed notebook retrieved **43 processed supplement
inputs** from verified article versions, including the named isoform container.
Instrument RAW and figure-data archives were excluded from this bounded review.
All original files, per-file URL/SHA256 receipts and published names remain in
the external archive under `source_review_inputs/`.

These are inputs for source/transform review, not automatic bindings to differently
named historical TSVs. Two sources remain version/media gaps. No arbitrary
supplement is renamed to pass the legacy filename gate. `review_inputs.json`
retains each file's original name, source PMCID, legacy requested filename,
license, URL, size, hash and association status. No installed values were changed.
