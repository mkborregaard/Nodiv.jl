# Standardized effect size (SOS metric) per grid cell from a simulation matrix.
function calculate_SOS(sims)
    sd = std.(eachcol(sims))
    me = mean.(eachcol(sims))
    (sims[1, :] .- me) ./ sd
end

# Default occupancy mask when the parent's occupancy is not supplied: cells whose
# descendant-richness column varies across the draws. The analysis entry points pass
# the deterministic, nsims-independent "parent clade present" mask instead.
_occupied(sims) = [!all(==(first(c)), c) for c in eachcol(sims)]

# Geographic node divergence (GND metric) from a simulation matrix, averaged over the
# OCCUPIED sites of the focal (parent) clade (Borregaard et al. 2014). Cells where the
# parent is absent carry no divergence signal and are excluded; a constant occupied
# column (descendant richness fixed across draws) contributes P ~ 1, i.e. no
# divergence. `occupied` is a per-cell boolean mask (default: the non-constant columns).
function calculate_GND(sims, occupied = _occupied(sims))
  # two internal convenience functions
  logit(p) = log(p / (1 - p))
  invlogit(p) = exp(p) / (1 + exp(p))

  n = size(sims, 1)
  idx = findall(occupied)
  isempty(idx) && return NaN
  r = [tiedrank(view(sims, :, j))[1] / (n + 1) for j in idx]
  # two-sided P (eqn 3); the -1/n keeps P off the 0/1 boundary so logit stays finite
  p = 1 .- 2 .* abs.(r .- 0.5) .- 1 / n
  α = mean(logit.(p))
  1 - invlogit(α)
end

# ---- effect-size alternatives to GND -----------------------------------------
# GND (eqn 4) combines per-cell two-sided P values in logit space. That makes it a
# p-value summary, not an effect size: its scale rides on `nsims` (P is bounded by
# ~1/nsims, so the maximum GND is 1 - O(1/nsims)), its no-divergence baseline is ~0.5
# rather than 0, and it saturates - the most divergent nodes pile against the ceiling.
# The functions below summarise the SOS field directly instead, giving a
# replication-stable magnitude in units of null SDs.

# Total spatial-divergence intensity: the root-mean-square SOS over the focal clade's
# occupied cells. Each SOS is ~standardised, so this is ~1 under the null and >1 under
# divergence. Replication-stable and not dominated by boundary cells. Constant occupied
# cells (sd = 0) contribute 0. Under the swap null, which fixes range size, a uniform
# richness/occupancy asymmetry between the clades is not flagged (use
# `calculate_GND_spatial` if you ever need to strip an offset explicitly).
function calculate_GND_rms(sims, occupied = _occupied(sims))
    sd = std.(eachcol(sims)); me = mean.(eachcol(sims))
    sos = (sims[1, :] .- me) ./ sd
    idx = findall(occupied)
    isempty(idx) && return NaN
    sqrt(mean(j -> isfinite(sos[j]) ? sos[j]^2 : 0.0, idx))
end

# Spatial-only intensity: as `calculate_GND_rms` but with the uniform offset
# (mean SOS = the clades' richness/occupancy asymmetry) removed, so it reflects only
# how over/under-representation varies ACROSS cells. Equals the SD of the SOS field.
function calculate_GND_spatial(sims, occupied = _occupied(sims))
    sd = std.(eachcol(sims)); me = mean.(eachcol(sims))
    sos = (sims[1, :] .- me) ./ sd
    idx = findall(occupied)
    isempty(idx) && return NaN
    std([isfinite(sos[j]) ? sos[j] : 0.0 for j in idx])
end

# Mean-square-SOS statistic of a richness row, standardised against the per-cell null
# moments, over occupied cells that actually vary. Internal helper for the SES and the
# Monte-Carlo P value, which share the null distribution of this statistic.
function _msos_null(sims, occupied)
    sd = std.(eachcol(sims)); me = mean.(eachcol(sims))
    keep = occupied .& (sd .> 0)
    any(keep) || return (NaN, Float64[])
    mek = me[keep]; sdk = sd[keep]
    msos(row) = (z = (view(row, keep) .- mek) ./ sdk; mean(abs2, z))
    msos(view(sims, 1, :)), [msos(view(sims, i, :)) for i in 2:size(sims, 1)]
end

# Standardised effect size of the divergence: the mean-square-SOS statistic of the
# empirical row, expressed in SDs of its null distribution. 0 = no divergence;
# replication-stable. NB the SES magnitude inflates with clade size (more cells ->
# tighter null), so prefer RMS-SOS for cross-node / cross-grain comparison.
function calculate_GND_ses(sims, occupied = _occupied(sims))
    Temp, Tnull = _msos_null(sims, occupied)
    isempty(Tnull) && return NaN
    s = std(Tnull)
    s == 0 ? NaN : (Temp - mean(Tnull)) / s
end

# Null-calibrated significance: the one-sided Monte-Carlo P value of the mean-square-
# SOS statistic (add-one estimator, bounded by 1/(nsims+1)). Small = divergent. NB this
# is a significance, so it carries a clade-size/power bias - larger clades clear a fixed
# P threshold more easily; `calculate_GND_rms` does not.
function calculate_GND_pval(sims, occupied = _occupied(sims))
    Temp, Tnull = _msos_null(sims, occupied)
    isempty(Tnull) && return NaN
    (1 + count(>=(Temp), Tnull)) / (length(Tnull) + 1)
end

# Same summaries from an already-computed per-cell SOS vector (e.g. a cached
# `NodeAnalysis.sos[node]`), so they can be derived without re-running the null.
# NaN/Inf cells (clade absent / constant) are dropped. `gnd_rms` = total intensity,
# `gnd_spatial` = spatial-only. The SES/P value need the null draws, so have no SOS-only form.
gnd_rms(sos::AbstractVector)     = (f = filter(isfinite, sos); isempty(f) ? NaN : sqrt(mean(abs2, f)))
gnd_spatial(sos::AbstractVector) = (f = filter(isfinite, sos); isempty(f) ? NaN : std(f))
