@testset "node analysis on toy data" begin
    assemblage, tree = toy_data()
    internal = [getnodename(tree, x) for x in traversal(tree, preorder) if !isleaf(tree, x)]
    others = setdiff(internal, TOY_ANALYSABLE)

    res = node_metrics(assemblage, tree; nsims=99)
    @test Set(keys(res.gnd)) == Set(internal)
    @test Set(keys(res.sos)) == Set(TOY_ANALYSABLE)
    for field in (:rms, :spatial, :ses, :pval, :varying)
        @test Set(keys(getfield(res, field))) == Set(TOY_ANALYSABLE)
    end
    @test all(isnan, res.gnd[n] for n in others)
    for n in TOY_ANALYSABLE
        @test length(res.sos[n]) == nsites(assemblage)
        @test isnan(res.sos[n][1])  # every species present: the null cannot vary
        @test isnan(res.sos[n][30])  # empty site
        @test all(isfinite, (res.gnd[n], res.rms[n], res.spatial[n], res.pval[n]))
    end
    # clades X and Y live on opposite halves of the grid
    @test res.rms["root"] > 1.5
    @test res.pval["root"] < 0.05
    @test res.nodes == internal
    @test sprint(show, res) == "NodeMetrics(11 internal nodes, 3 analysed, 30 cells)"
    for n in TOY_ANALYSABLE
        # rms and spatial summarise the cells where the null varies, so they follow
        # from the SOS alone
        @test res.rms[n] == sos_rms(res.sos[n])
        @test res.spatial[n] == sos_sd(res.sos[n])
        @test 0 < res.varying[n] < 1  # site 1 is occupied but cannot vary
    end

    ana = node_analysis(assemblage, tree; nsims=99)
    @test ana.nodes == internal
    @test Set(keys(ana.gnd)) == Set(internal)
    @test Set(keys(ana.sos)) == Set(TOY_ANALYSABLE)
    @test isnan.(ana.sos["X"]) == isnan.(res.sos["X"])

    sos_scores, gnd = process_node(assemblage, tree, "root"; nsims=99)
    @test length(sos_scores) == nsites(assemblage)
    @test isfinite(gnd)
    sos_scores, gnd = process_node(assemblage, tree, first(others); nsims=99)
    @test all(isnan, sos_scores)
    @test isnan(gnd)

    sims = simulate_descendants(get_clade(assemblage, tree, "root"), tree, "X"; nsims=20)
    @test size(sims) == (21, nsites(assemblage))
    @test sims[1, :] == richness(get_clade(assemblage, tree, "X"))
end

@testset "the varying share links the two readings of RMS-SOS" begin
    assemblage, tree = toy_data()
    clade = get_clade(assemblage, tree, "root")
    occupied = richness(clade) .> 0
    sims = simulate_descendants(clade, tree, "X"; nsims=99)
    sos_scores = sos(sims)
    share = Nodiv._varying_share(sos_scores, occupied)
    @test share ≈ count(isfinite, sos_scores) / count(occupied)
    @test share < 1
    # RMS-SOS over all occupied cells, the fixed ones counted as 0
    squares = [isfinite(s) ? s^2 : 0.0 for s in sos_scores[occupied]]
    @test sqrt(mean(squares)) ≈ sqrt(share) * sos_rms(sos_scores)
end

@testset "divergent_nodes" begin
    assemblage, tree = toy_data()
    res = node_metrics(assemblage, tree; nsims=99)
    above(d, t) = Set(n for (n, v) in d if !isnan(v) && v > t)

    @test Set(divergent_nodes(res)) == above(res.rms, 1.5)
    @test Set(divergent_nodes(res; threshold=0)) == Set(TOY_ANALYSABLE)
    @test Set(divergent_nodes(res; by=:gnd)) == above(res.gnd, 0.8)
    @test Set(divergent_nodes(res; by=:pval)) == Set(n for (n, p) in res.pval if p < 0.05)
    @test Set(divergent_nodes(res.gnd; threshold=0.5)) == above(res.gnd, 0.5)
    @test_throws ArgumentError divergent_nodes(res; by=:ses)

    # the result types give the nodes in tree order
    divergent = divergent_nodes(res; threshold=0)
    @test divergent == filter(in(Set(TOY_ANALYSABLE)), res.nodes)
    @test divergent_nodes(res; by=:pval, threshold=1.1) == divergent
    # a plain Dict has no tree: most divergent first
    ranked = divergent_nodes(res.rms; threshold=0)
    @test issorted([res.rms[n] for n in ranked]; rev=true)

    ana = node_analysis(assemblage, tree; nsims=99)
    @test Set(divergent_nodes(ana; threshold=0.5)) == above(ana.gnd, 0.5)
    @test divergent_nodes(ana; threshold=0) == filter(in(Set(TOY_ANALYSABLE)), ana.nodes)
end

@testset "prune_to_shared!" begin
    assemblage, _ = toy_data()
    tree = parsenewick("(" * chop(TOY_NEWICK) * ":1,m:1)top;")
    @test "m" in getleafnames(tree)
    @test prune_to_shared!(tree, assemblage) === tree
    @test sort(getleafnames(tree)) == string.('a':'l')
end
