# CLAUDE.md — Nodiv

## What this is

Nodiv is a local development Julia package (`]dev`), authored by Michael Borregaard and co-developed with Claude from earlier R work. It implements the node-based SOS/GND analysis of Borregaard et al. 2014 (the paper PDF is available in the project). Edits here are real changes to a real package, not throwaway exploration.

Package relationships:
- Nodiv depends on **SpatialEcology.jl** (also authored here) for the `Assemblage` type and spatial plumbing.
- Nodiv is driven and exercised by **NodivWorkshop** (a separate directory/repo) whose `script.jl` is the integration workflow. Nodiv is the package; the workshop is the application that uses it.
- Edit Nodiv source freely. Treat SpatialEcology and NodivWorkshop as separate boundaries — call out explicitly when a change would require touching either of them.

## Running and testing

*(Verify these against the actual repo before relying on them — inferred from standard Julia layout.)*
- Work in the package env: `julia --project=.`
- Run the test suite: `julia --project=. -e 'using Pkg; Pkg.test()'` (confirm `test/runtests.jl` exists and passes first).

## Expensive step — do not recompute

The integration check is `NodivWorkshop/script.jl`. It contains one very expensive step: `node_metrics(assemblage, tree; nsims=...)` runs the null-model randomizations over the whole tree (~18k cells for the geographic scan) for both spaces. Its results are cached to the workshop's `data/node_analysis.jld2` as `res_e` and `res_g`.

Do **not** trigger a recompute. Work against the cached `res_e` / `res_g`, which load from that file. If the cache is absent or you think it needs regenerating, ask first — it is a long run.

## Constraints

- Keep the public signature of `sos_distances(res, nodes; …)` stable. `NodivWorkshop/script.jl` calls it and must keep working. Add keyword arguments with defaults; do not change the positional arguments or the return type (a pairwise distance matrix consumed by `fit(MDS, …; distances = true)`).
- This is a published method. Changes to metric or distance definitions must be scientifically deliberate and match the design spec, not convenience refactors.

## Current task

See `docs/sos_pattern_grouping_design.md` for the active design spec: grouping divergent nodes by SOS-pattern similarity. That document is the source of truth for what to build and why. Start there.
