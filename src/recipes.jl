# Plot a node in a 2x2 grid: the parent clade (top-left), the descendant SOS
# mapped over the assemblage's coordinates (top-right, RdYlBu), and the two child
# clades on the bottom row. Defined as a plot recipe so the package depends only
# on RecipesBase; the plotting backend (Plots) is supplied by the caller. Use as
# `plot_node(assemblage, tree, node)`.
@userplot Plot_Node

@recipe function f(pn::Plot_Node)
    assemblage, tree, node = pn.args[1:3]
    ch1, ch2 = getchildren(tree, node)[1:2]
    assm = get_clade(assemblage, tree, node)
    assmch1 = get_clade(assm, tree, ch1)
    assmch2 = get_clade(assm, tree, ch2)

    # SOS for the top-right panel, taken from the cached analysis result supplied as
    # the 4th argument - either a `NodeAnalysis`/`NodeMetrics` (looked up by node) or a
    # precomputed SOS vector. Pass the result of `node_metrics`/`node_analysis`; the
    # panel is never recomputed on the fly.
    length(pn.args) >= 4 ||
        error("plot_node needs the analysis result (or a precomputed SOS vector) as the " *
              "4th argument, e.g. plot_node(assemblage, tree, node, res)")
    cached = pn.args[4]
    sos = cached isa Union{NodeAnalysis, NodeMetrics} ? cached.sos[node] : cached

    layout := (2, 2)
    size --> (900, 800)

    @series begin              # top-left: parent clade
        subplot := 1
        title := string(node)
        assm
    end
    @series begin              # top-right: SOS in environmental space
        subplot := 2
        title := "SOS"
        fillcolor := :RdYlBu
        clim := (-8, 8)
        sos, assm
    end
    @series begin              # bottom-left: child 1
        subplot := 3
        title := string(ch1)
        assmch1
    end
    @series begin              # bottom-right: child 2
        subplot := 4
        title := string(ch2)
        assmch2
    end
end

# Plot per-node values (e.g. GND) on the tree. `gndvals` is a Dict of node name
# => value; a coloured marker is drawn at every node present in it (NaN values
# skipped) and nothing elsewhere. Pass a NodeAnalysis (or its `gnd` Dict) to show
# all analysable nodes, or a filtered Dict (e.g. only divergent nodes) for a subset.
# Use as `plot_gnd(tree, gnd)`.
@userplot Plot_Gnd

@recipe function f(pg::Plot_Gnd)
    tree, gndvals = pg.args
    gndvals isa Union{NodeAnalysis, NodeMetrics} && (gndvals = gndvals.gnd)
    shown = Dict(k => v for (k, v) in gndvals if !isnan(v))

    # The tree recipe draws a marker at every node and renders NaN-marker_z nodes
    # solid, so hide the non-shown nodes with size 0. markersize must follow the
    # recipe's node order, which is Phylo's layout helper _findxy.
    layoutnodes = Phylo._findxy(tree)[3]

    treetype --> :fan
    showtips --> false
    markerstrokewidth --> 0
    color --> :YlOrRd
    clim --> (0, 1)
    size --> (1000, 1000)
    # `markersize` sets the size of the *shown* markers (default 6) and is overridable -
    # pass a scalar; it is expanded here so the hidden (non-shown) nodes always stay at 0.
    markersize --> 6
    ms = plotattributes[:markersize]
    marker_z := shown
    markersize := [haskey(shown, node) ? ms : 0 for node in layoutnodes]
    tree
end
