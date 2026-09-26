# Per-node divergence scores from a simulation matrix `sims` (see `simulate_descendants`):
# row 1 holds the observed richness of the descendant clade in each cell, and each further
# row one draw of the null model.

# Per-cell mean and standard deviation of descendant richness over all rows of `sims`
# (the observed row and the null draws), computed once and shared by the metrics.
function _moments(sims)
    me = mean.(eachcol(sims))
    sd = [std(column; mean=m) for (column, m) in zip(eachcol(sims), me)]
    return me, sd
end

_sos(sims, me, sd) = (view(sims, 1, :) .- me) ./ sd

"""
    sos(sims) -> Vector{Float64}

The specific overrepresentation score (SOS) of every cell: the observed richness of the
descendant clade, standardised against the null model as `(observed - mean) / sd` over
the rows of `sims`. Positive where the descendant clade is over-represented.

`sims` is the matrix from [`simulate_descendants`](@ref): the observed richness in row 1
and one null draw per further row, one column per cell. Cells where the richness does
not vary have no SOS (`NaN`): cells where the parent clade is absent, and occupied cells
where the null model cannot vary, e.g. because every species of the parent clade is there.
"""
sos(sims::AbstractMatrix) = _sos(sims, _moments(sims)...)

# Default occupancy mask when the parent's occupancy is not supplied: cells whose
# descendant-richness column varies across the draws. The analysis entry points pass
# the deterministic, nsims-independent "parent clade present" mask instead.
_occupied(sims) = [!all(==(first(c)), c) for c in eachcol(sims)]

"""
    gnd(sims[, occupied]) -> Float64

The geographic node divergence (GND) of Borregaard et al. (2014, eqns 3-4). Each cell's
two-sided P value of the observed richness among the null draws is combined in logit space
over the cells in `occupied`, as `1 - invlogit(mean(logit(P)))`. `occupied` is a mask
with one element per cell, normally where the parent clade is present; by default the
cells whose richness varies. Occupied cells where the null model cannot vary get P ≈ 1,
so count as no divergence. `NaN` if no cell is occupied.

GND summarises P values rather than being an effect size: its no-divergence baseline is
about 0.5, its maximum is `1 - O(1/nsims)`, and the most divergent nodes pile up near that
ceiling. [`sos_rms`](@ref) is the recommended score.
"""
function gnd(sims::AbstractMatrix, occupied=_occupied(sims))
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

"""
    sos_rms(sims) -> Float64
    sos_rms(sos_scores::AbstractVector) -> Float64

RMS-SOS, the recommended divergence score: the root-mean-square [`sos`](@ref) over the
cells where the null model varies, from a simulation matrix or from SOS values (e.g. a
cached `res.sos[node]`). As each SOS is standardised, it is about 1 without divergence and
larger with it, in units of null standard deviations, and it does not depend on `nsims`.
Cells without an SOS (`NaN`) are left out; `NaN` if there are none.

The cells left out include occupied cells where the null model cannot vary, such as
cells holding every species of the parent clade. Where such cells are much of a clade's
range, the score rests on the rest of it; `node_metrics` reports their share as `varying`.
"""
sos_rms(sims::AbstractMatrix) = sos_rms(sos(sims))
function sos_rms(sos_scores::AbstractVector)
    finite = filter(isfinite, sos_scores)
    return isempty(finite) ? NaN : sqrt(mean(abs2, finite))
end

"""
    sos_sd(sims) -> Float64
    sos_sd(sos_scores::AbstractVector) -> Float64

The standard deviation of the [`sos`](@ref) over the cells where the null model varies.
Like [`sos_rms`](@ref), but without any uniform offset between the two descendant clades,
so it measures only how their relative representation changes across cells.
"""
sos_sd(sims::AbstractMatrix) = sos_sd(sos(sims))
function sos_sd(sos_scores::AbstractVector)
    finite = filter(isfinite, sos_scores)
    return isempty(finite) ? NaN : std(finite)
end

# Mean-square-SOS statistic of each richness row, standardised against the per-cell
# moments, over the cells in `keep`: the observed row's value and the null rows' values.
# Internal helper for the SES and the Monte-Carlo P value, which share the null
# distribution of this statistic.
function _msos_null(sims, me, sd, keep)
    any(keep) || return (NaN, Float64[])
    mek = me[keep]
    sdk = sd[keep]
    msos(row) = mean(abs2, (view(row, keep) .- mek) ./ sdk)
    return msos(view(sims, 1, :)), [msos(view(sims, i, :)) for i in 2:size(sims, 1)]
end

# Over the cells in `occupied` that vary
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

"""
    divergence_ses(sims) -> Float64

The standardised effect size of the divergence: the mean square [`sos`](@ref) of the
observed row, over the cells where the null model varies, in standard deviations of its
distribution among the null draws. About 0 without divergence.

Its magnitude grows with clade size, as more cells give a tighter null distribution, so do
not compare it across nodes or grains; use [`sos_rms`](@ref) for that.
"""
divergence_ses(sims::AbstractMatrix) = _ses(_msos_null(sims, trues(size(sims, 2)))...)

"""
    divergence_pval(sims) -> Float64

The Monte Carlo P value of the divergence: the share of null draws whose mean square
[`sos`](@ref) (over the cells where the null model varies) is at least that of the
observed row, by the add-one estimator, so at least `1 / (nsims + 1)`. Small means
divergent. Being a significance, it depends on power: larger clades clear a fixed
threshold more easily.
"""
divergence_pval(sims::AbstractMatrix) = _pval(_msos_null(sims, trues(size(sims, 2)))...)

# Share of the focal clade's occupied cells where the null model varies, i.e. where the
# SOS is defined. The RMS-SOS over all occupied cells, counting the others as 0, is
# `sqrt(share) * sos_rms(sos_scores)`.
function _varying_share(sos_scores, occupied)
    n = count(occupied)
    n == 0 && return NaN
    return count(i -> occupied[i] && isfinite(sos_scores[i]), eachindex(sos_scores)) / n
end
