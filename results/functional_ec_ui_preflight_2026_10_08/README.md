# Functional UI packaging preflight

The first wrapper failed before writing any output: `ValueError` from applying
`Path.relative_to()` to a relative pilot path and an absolute repository root.
The source artifact was unchanged. Resolve the pilot path before recording the
relative source location, then execute a new packaging output. Future wrapper
failures preserve an executed diagnostic notebook rather than relying on this
record of the first preflight error.
