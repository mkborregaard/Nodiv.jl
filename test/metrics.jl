@testset "divergence metrics from a sims matrix" begin
    Random.seed!(42)
    # 100 null draws + an empirical first row; no divergence
    null = Float64.(rand(0:4, 101, 40))
    # same, but the empirical row is over-represented in 12 cells -> divergence
    div = copy(null)
    div[1, 1:12] .= 12.0

    # RMS-SOS is ~1 under the null and clearly >1 under divergence
    @test isapprox(calculate_GND_rms(null), 1.0; atol=0.4)
    @test calculate_GND_rms(div) > 2

    # GND baseline is ~0.5 (not 0) under the null and rises under divergence
    @test 0.3 < calculate_GND(null) < 0.7
    @test calculate_GND(div) > calculate_GND(null)

    # Monte-Carlo P value: in [1/(n+1), 1], smaller for the divergent case
    @test 1 / 101 <= calculate_GND_pval(div) <= 1
    @test calculate_GND_pval(div) <= calculate_GND_pval(null)

    # SES rises with divergence
    @test calculate_GND_ses(div) > calculate_GND_ses(null)
end

@testset "occupied mask is deterministic and excludes empty cells" begin
    Random.seed!(7)
    sims = Float64.(rand(0:4, 51, 8))
    sims[1, 1:3] .= 10.0
    empty = hcat(sims, zeros(51))  # a parent-absent (all-zero) cell
    occ = vcat(trues(8), false)
    # adding an unoccupied cell must not change the score when masked out
    @test isapprox(calculate_GND_rms(empty, occ), calculate_GND_rms(sims); rtol=1e-9)
    @test isapprox(calculate_GND(empty, occ), calculate_GND(sims); rtol=1e-9)
    # a constant occupied column contributes 0 to RMS, not NaN
    flat = hcat(sims, fill(2.0, 51))
    @test isfinite(calculate_GND_rms(flat, vcat(trues(8), true)))
end

@testset "cached-SOS helpers match the sims-level functions" begin
    Random.seed!(1)
    sims = Float64.(rand(0:5, 81, 30))
    sims[1, 1:6] .= 9.0
    sos_scores = calculate_SOS(sims)
    @test isapprox(gnd_rms(sos_scores), calculate_GND_rms(sims); rtol=1e-9)
    @test isapprox(gnd_spatial(sos_scores), calculate_GND_spatial(sims); rtol=1e-9)
end
