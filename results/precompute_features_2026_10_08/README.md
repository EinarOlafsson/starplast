# Frozen numeric feature operator

Executed `pilot.ipynb` verifies a reusable numeric operator against the original
FN-EC-01 artifact. All 1,226 ordered entities, 679 training entities, 349 permitted
columns and missing-value masks are retained. Training ranks, frozen query ECDFs
and historical shared training vectors/distributions agree exactly.

The typed child binds the operator identity. A deliberate one-job stop preserves
the first result; resume builds only the child; replay uses no builders. Source,
split, exclusions, settings and guard-code receipts accompany the outputs.
`manifest.json` records the 33 executed output files; this explanatory README was
added after execution and is outside that receipt.

The run took 16.64 seconds and recorded 19,999,824 evidence bytes before this README.
The external cgroup limit was 400 MiB; observed process VmHWM was 401.6 MiB
(421,117,952 bytes). Process RSS and cgroup accounting differ; this is not evidence
of a process peak below 400 MiB. Cooperative checkpoints, 120-second/two-attempt
limits and 32 MiB artifact/package limits are recorded in the runner reports.

This verifies software reuse, not biological accuracy. No new classifier,
calibrated confidence or deployment output was produced. Graph, representation,
fitting and deployment builders and broader organism coverage remain unfinished.
