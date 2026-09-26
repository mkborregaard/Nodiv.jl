"""
    nodespecies(tree, node) -> Vector{String}

The names of the tips (species) descending from `node`, or the node itself if it is a tip.
"""
function nodespecies(tree, node)
    # A concrete Vector{String}: an empty or tip result typed like `getdescendants`
    # (Union{Nothing,String}) misses the name-matching method of `view`
    isleaf(tree, node) && return String[getnodename(tree, node)]
    return String[x for x in getdescendants(tree, node) if isleaf(tree, x)]
end

"""
    get_clade(assemblage, tree, node)

A view of `assemblage` with only the species descending from `node`.
"""
get_clade(assemblage, tree, node) = view(assemblage; species=nodespecies(tree, node))

"""
    prune_to_shared!(tree, assemblages...) -> tree

Remove the tips of `tree` that are missing from any of `assemblages`, so every species in
the tree has data. Returns the pruned tree.
"""
function prune_to_shared!(tree, assemblages...)
    shared = intersect(getleafnames(tree), speciesnames.(assemblages)...)
    keeptips!(tree, shared)
    return tree
end
