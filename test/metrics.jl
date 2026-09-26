@testset "divergence metrics from a sims matrix" begin
    Random.seed!(42)
    # 100 null draws + an observed first row; no divergence
    null = Float64.(rand(0:4, 101, 40))
    # same, but the observed row is over-represented in 12 cells -> divergence
    div = copy(null)
    div[1, 1:12] .= 12.0

    # RMS-SOS is ~1 under the null and clearly >1 under divergence
    @test isapprox(sos_rms(null), 1.0; atol=0.4)
    @test sos_rms(div) > 2

    # GND baseline is ~0.5 (not 0) under the null and rises under divergence
    @test 0.3 < gnd(null) < 0.7
    @test gnd(div) > gnd(null)

    # Monte-Carlo P value: in [1/(n+1), 1], smaller for the divergent case
    @test 1 / 101 <= divergence_pval(div) <= 1
    @test divergence_pval(div) <= divergence_pval(null)

    # SES rises with divergence
    @test divergence_ses(div) > divergence_ses(null)
end

@testset "cells without an SOS are left out" begin
    Random.seed!(7)
    sims = Float64.(rand(0:4, 51, 8))
    sims[1, 1:3] .= 10.0
    # an empty cell (parent absent) and a fixed occupied cell (the null cannot vary)
    padded = hcat(sims, zeros(51), fill(2.0, 51))
    @test all(isnan, sos(padded)[9:10])
    for score in (sos_rms, sos_sd, divergence_ses, divergence_pval)
        @test isapprox(score(padded), score(sims); rtol=1e-9)
    end
    # GND ignores unoccupied cells and counts a fixed occupied cell as no divergence
    @test isapprox(gnd(padded, vcat(trues(8), false, false)), gnd(sims); rtol=1e-9)
    @test gnd(padded, vcat(trues(8), false, true)) < gnd(sims)
end

@testset "SOS values and sims give the same summaries" begin
    Random.seed!(1)
    sims = Float64.(rand(0:5, 81, 30))
    sims[1, 1:6] .= 9.0
    @test sos_rms(sos(sims)) == sos_rms(sims)
    @test sos_sd(sos(sims)) == sos_sd(sims)
    @test isnan(sos_rms([NaN, NaN]))
end
