# CLAUDE.md — Nodiv

## What this is

Nodiv is a local development Julia package (`]dev`), authored by Michael Borregaard and co-developed with Claude from earlier R work. It implements the node-based SOS/GND analysis of Borregaard et al. 2014 (the paper PDF is available in the project). Edits here are real changes to a real package, not throwaway exploration.

Package relationships:
- Nodiv depends on **SpatialEcology.jl** (also authored here) for the `Assemblage` type and spatial plumbing.
- Nodiv is driven and exercised by **NodivWorkshop** (a separate directory/repo) whose `script.jl` is the integration workflow. Nodiv is the package; the workshop is the application that uses it.
- Edit Nodiv source freely. Treat SpatialEcology and NodivWorkshop as separate boundaries — call out explicitly when a change would require touching either of them.

## Running and testing

- Test deps are in `[extras]`/`[targets]` of `Project.toml` (no `test/Project.toml`); `julia_eval` with `env_path = <root>/test/` activates them, then `include("test/runtests.jl")`.
- Full run as on CI, with deprecation warnings checked: `julia --project=. -e 'using Pkg; Pkg.test(; julia_args=["--depwarn=yes"])'`.
- Tests run on a toy assemblage (`test/toydata.jl`) in seconds; they never need the workshop data.
- Code is BlueStyle, enforced by JuliaFormatter 2.14 (`.JuliaFormatter.toml`, checked in CI by `Format.yml`). Run `JuliaFormatter.format(".")` before committing.

## Expensive step — do not recompute

The integration check is `NodivWorkshop/script.jl`. It contains one very expensive step: `node_metrics(assemblage, tree; nsims=...)` runs the null-model randomizations over the whole tree (~18k cells for the geographic scan) for both spaces. Its results are cached to the workshop's `data/node_analysis.jld2` as `res_e` and `res_g`.

Do **not** trigger a recompute casually. Work against the cached `res_e` / `res_g`, which load from that file. The cache is regenerated whenever a change actually invalidates it (e.g. a change to the result types or to how a stored score is defined); say so when a change does, and ask before starting the long run yourself.

## Constraints

- Keep the public signature of `sos_distances(res, nodes; …)` stable. `NodivWorkshop/script.jl` calls it and must keep working. Add keyword arguments with defaults; do not change the positional arguments or the return type (a pairwise distance matrix consumed by `sos_ordination`, `sos_clusters` and `sos_similarity_communities`).
- The Plots recipes (`plot_gnd`, `plot_node`) stay in Nodiv for users who do not use Makie. Do not move or deprecate them unless the user explicitly says so.
- New function names use underscores for multi-word names and are never capitalized (types keep CamelCase).
- Division of labour: general analysis infrastructure goes in Nodiv, visualisation in NodivMakie, anything specific to the workshop's data (birds, AVONET traits, genus/family/order names) in NodivWorkshop. Nodiv's docs must not point to the workshop.
- This is a published method. Changes to metric or distance definitions must be scientifically deliberate and match the design spec, not convenience refactors.

## Design spec

The design spec for grouping divergent nodes by SOS-pattern similarity (`sos_distances`, `sos_ordination`, `sos_clusters`, `sos_similarity_communities`) is `docs/sos_pattern_grouping_design.md` in NodivWorkshop, as it is written about the workshop's data. It is the source of truth for what these functions do and why.
