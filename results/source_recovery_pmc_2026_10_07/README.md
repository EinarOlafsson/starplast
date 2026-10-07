# Initial cloud metadata diagnostic — superseded

The first cloud probe retained primary version listings and JSON responses.
Its conservative parser did not yet accept boolean manuscript/retraction flags
and S3 URLs. It recovered no inputs. Original parser code is under `code/`;
the manifest and executed notebook retain this failure without fabricating
missing-file evidence. The corrected tested run is the sibling `_v2` snapshot.
