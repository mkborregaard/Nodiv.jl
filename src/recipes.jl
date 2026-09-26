# Plots recipes. Nodiv depends only on RecipesBase; the caller loads Plots.

"""
    plot_node(assemblage, tree, node, res)

Plot a node with Plots, in a 2 x 2 grid: the richness of its clade (top left), its SOS
map (top right), and the richness of its two descendant clades (bottom row). `res` is the
result of [`node_metrics`](@ref) or [`node_analysis`](@ref), or the node's SOS vector;
nothing is recomputed.
"""
@userplot Plot_Node

@recipe function f(pn::Plot_Node)
    assemblage, tree, node = pn.args[1:3]
    ch1, ch2 = getchildren(tree, node)[1:2]
    assm = get_clade(assemblage, tree, node)
    assmch1 = get_clade(assm, tree, ch1)
    assmch2 = get_clade(assm, tree, ch2)

    if length(pn.args) < 4
        msg =
            "plot_node needs the analysis result (or a precomputed SOS vector) as " *
            "the 4th argument, e.g. plot_node(assemblage, tree, node, res)"
        throw(ArgumentError(msg))
    end
    cached = pn.args[4]
    sos_scores = cached isa AbstractNodeResult ? cached.sos[node] : cached

    layout := (2, 2)
    size --> (900, 800)

    @series begin              # top-left: parent clade
        subplot := 1
        title := string(node)
        assm
    end
    @series begin              # top-right: SOS
        subplot := 2
        title := "SOS"
        fillcolor := :RdYlBu
        clim := (-8, 8)
        sos_scores, assm
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

"""
    plot_gnd(tree, res)
    plot_gnd(tree, values::AbstractDict)

Plot per-node values on the tree with Plots, as coloured markers: by default the GND of
the analysed nodes of `res`. Pass a Dict of node name => value instead, e.g.
`res.rms` or only the divergent nodes, to show other values or a subset. Nodes with no
value or a `NaN` get no marker. Draws a fan tree coloured `:YlOrRd` over `(0, 1)`;
`markersize` sets the size of the shown markers.
"""
@userplot Plot_Gnd

@recipe function f(pg::Plot_Gnd)
    tree, gndvals = pg.args
    gndvals isa AbstractNodeResult && (gndvals = gndvals.gnd)
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
