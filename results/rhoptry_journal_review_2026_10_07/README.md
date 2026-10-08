# Initial journal retrieval and workbook diagnostic — superseded

The public publisher's Dataset EV1 archive downloaded and passed ZIP integrity
checks. The initial parser then tried to read a macOS resource-fork entry named
like an XLSX workbook. That is not an actual workbook, and the failure is retained
in `review.json`, the executed notebook and original code copy.

The corrected `_v2` review excludes resource forks and verifies/reuses this exact
archive through the original receipt. No file is redownloaded or overwritten to
hide the failure. This directory alone does not establish source equivalence.
