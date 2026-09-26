# Per-cell mean and standard deviation of descendant richness over all rows of `sims`
# (the empirical row and the null draws), computed once and shared by the metrics.
function _moments(sims)
    me = mean.(eachcol(sims))
    sd = [std(column; mean=m) for (column, m) in zip(eachcol(sims), me)]
    return me, sd
end

_sos(sims, me, sd) = (view(sims, 1, :) .- me) ./ sd

# Standardized effect size (SOS metric) per grid cell from a simulation matrix.
calculate_SOS(sims) = _sos(sims, _moments(sims)...)

# Default occupancy mask when the parent's occupancy is not supplied: cells whose
# descendant-richness column varies across the draws. The analysis entry points pass
# the deterministic, nsims-independent "parent clade present" mask instead.
_occupied(sims) = [!all(==(first(c)), c) for c in eachcol(sims)]

# Geographic node divergence (GND metric) from a simulation matrix, averaged over the
# OCCUPIED sites of the focal (parent) clade (Borregaard et al. 2014). Cells where the
# parent is absent carry no divergence signal and are excluded; a constant occupied
# column (descendant richness fixed across draws) contributes P ~ 1, i.e. no
# divergence. `occupied` is a per-cell boolean mask (default: the non-constant columns).
function calculate_GND(sims, occupied=_occupied(sims))
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
    return 1 - invlogit(α)
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
function calculate_GND_rms(sims, occupied=_occupied(sims))
    sos = calculate_SOS(sims)
    idx = findall(occupied)
    isempty(idx) && return NaN
    return sqrt(mean(j -> isfinite(sos[j]) ? sos[j]^2 : 0.0, idx))
end

# Spatial-only intensity: as `calculate_GND_rms` but with the uniform offset
# (mean SOS = the clades' richness/occupancy asymmetry) removed, so it reflects only
# how over/under-representation varies ACROSS cells. Equals the SD of the SOS field.
function calculate_GND_spatial(sims, occupied=_occupied(sims))
    sos = calculate_SOS(sims)
    idx = findall(occupied)
    isempty(idx) && return NaN
    return std([isfinite(sos[j]) ? sos[j] : 0.0 for j in idx])
end

# Mean-square-SOS statistic of each richness row, standardised against the per-cell
# moments, over the cells in `keep`: the empirical row's value and the null rows' values.
# Internal helper for the SES and the Monte-Carlo P value, which share the null
# distribution of this statistic.
function _msos_null(sims, me, sd, keep)
    any(keep) || return (NaN, Float64[])
    mek = me[keep]
    sdk = sd[keep]
    msos(row) = mean(abs2, (view(row, keep) .- mek) ./ sdk)
    return msos(view(sims, 1, :)), [msos(view(sims, i, :)) for i in 2:size(sims, 1)]
end

# Over the occupied cells that actually vary
function _msos_null(sims, occupied)
    me, sd = _moments(sims)
    return _msos_null(sims, me, sd, occupied .& (sd .> 0))
end

function _ses(stat, nullstats)
    isempty(nullstats) && return NaN
    s = std(nullstats)
    return s == 0 ? NaN : (stat - mean(nullstats)) / s
end

function _pval(stat, nullstats)
    isempty(nullstats) && return NaN
    return (1 + count(>=(stat), nullstats)) / (length(nullstats) + 1)
end

# Standardised effect size of the divergence: the mean-square-SOS statistic of the
# empirical row, expressed in SDs of its null distribution. 0 = no divergence;
# replication-stable. NB the SES magnitude inflates with clade size (more cells ->
# tighter null), so prefer RMS-SOS for cross-node / cross-grain comparison.
calculate_GND_ses(sims, occupied=_occupied(sims)) = _ses(_msos_null(sims, occupied)...)

# Null-calibrated significance: the one-sided Monte-Carlo P value of the mean-square-
# SOS statistic (add-one estimator, bounded by 1/(nsims+1)). Small = divergent. NB this
# is a significance, so it carries a clade-size/power bias - larger clades clear a fixed
# P threshold more easily; `calculate_GND_rms` does not.
calculate_GND_pval(sims, occupied=_occupied(sims)) = _pval(_msos_null(sims, occupied)...)

# Same summaries from an already-computed per-cell SOS vector (e.g. a cached
# `NodeAnalysis.sos[node]`), so they can be derived without re-running the null.
# NaN/Inf cells (clade absent / constant) are dropped. `gnd_rms` = total intensity,
# `gnd_spatial` = spatial-only. The SES/P value need the null draws, so have no SOS-only form.
function gnd_rms(sos::AbstractVector)
    return (f=filter(isfinite, sos); isempty(f) ? NaN : sqrt(mean(abs2, f)))
end
gnd_spatial(sos::AbstractVector) = (f=filter(isfinite, sos); isempty(f) ? NaN : std(f))

# Share of the focal clade's occupied cells where the null model varies, i.e. where the
# SOS is defined. The RMS-SOS over all occupied cells, counting the others as 0, is
# `sqrt(share) * gnd_rms(sos)`.
function _varying_share(sos, occupied)
    n = count(occupied)
    n == 0 && return NaN
    return count(i -> occupied[i] && isfinite(sos[i]), eachindex(sos)) / n
end
