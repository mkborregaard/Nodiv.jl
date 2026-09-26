# All tips/species descending from a node (or the node itself if it is a tip).
# Returns a concretely-typed Vector{String}: `getdescendants` is typed as
# Union{Nothing,String}, and an empty/tip result of that type misses the name
# matching method used by the downstream `view`.
function nodespecies(tree, node)
    isleaf(tree, node) && return String[getnodename(tree, node)]
    return String[x for x in getdescendants(tree, node) if isleaf(tree, x)]
end

# Subset an Assemblage to the clade descending from a node.
get_clade(assemblage, tree, node) = view(assemblage; species=nodespecies(tree, node))

# Prune `tree` in place to the tips it shares with all the given assemblage(s),
# so clade subsetting never references a species absent from the data. Returns
# the tree.
function prune_to_shared!(tree, assemblages...)
    shared = intersect(getleafnames(tree), speciesnames.(assemblages)...)
    keeptips!(tree, shared)
    return tree
end
