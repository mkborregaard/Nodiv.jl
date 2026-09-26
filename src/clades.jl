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
    clade_richness(assemblage, tree) -> Function
    clade_richness(assemblage, tree, node) -> Vector{Int}

The richness of the clade of `node` in each cell of `assemblage`, the same as
`richness(get_clade(assemblage, tree, node))` but fast. Species of the clade missing from
`assemblage` are left out.

The two-argument form indexes the assemblage once and returns a function of the node
name, for the richness of many clades: on a large assemblage (10,000 species, 18,000 cells)
each clade then takes milliseconds instead of most of a second.
"""
function clade_richness(assemblage, tree)
    presence = sparse(permutedims(occurrences(assemblage) .> 0))  # Sites x species
    index = Dict(String(sp) => j for (j, sp) in enumerate(speciesnames(assemblage)))
    rows, vals = rowvals(presence), nonzeros(presence)
    return function (node)
        r = zeros(Int, size(presence, 1))
        for sp in nodespecies(tree, node)
            j = get(index, sp, 0)
            j == 0 && continue
            for k in nzrange(presence, j)
                if vals[k]
                    r[rows[k]] += 1
                end
            end
        end
        return r
    end
end

clade_richness(assemblage, tree, node) = clade_richness(assemblage, tree)(node)

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
