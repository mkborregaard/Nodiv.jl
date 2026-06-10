using CSV, DataFrames, SpatialEcology, Phylo, Plots
using MultivariateStats, Statistics

using Nodiv

const datadir = "/Users/cvg147/Dropbox/Arbejde/Current projects/NodivWorkshop"

### Load data and create objects

pam = CSV.read(joinpath(datadir, "PAM_E.csv"), DataFrame)

env = CSV.read(joinpath(datadir, "Env.csv"), DataFrame)

# Cross.csv links the PAM/BirdLife taxonomy (Species1) to the tree/BirdTree
# taxonomy (Species3). Use it to relabel PAM species to their tree-tip names so
# the datasets share as many species as possible; unmatched names fall back to a
# plain underscore conversion. (Several BirdLife taxa can map to one BirdTree
# tip, in which case their occurrences are merged into that tip.)
cross = CSV.read(joinpath(datadir, "Cross.csv"), DataFrame)
namemap = Dict(string(r.Species1) => replace(string(r.Species3), " " => "_")
               for r in eachrow(cross) if !ismissing(r.Species1) && !ismissing(r.Species3))

# PAM_E is long-format presence/absence: ID_env links to the env cell (site),
# Species is the species name, and there is no abundance column. Build the
# 3-column phylocom table [site, abundance, species] that ComMatrix expects.
phylocom = DataFrame(site = string.(pam.ID_env),
                     abundance = 1,
                     species = [get(namemap, s, replace(s, " " => "_")) for s in pam.Species])

# Coordinates live in environmental (PC) space: use each cell's bin midpoint,
# x = mean(xmin, xmax), y = mean(ymin, ymax). SpatialEcology aligns the matrix
# and coordinates by row order (not by name), so reduce + reorder env to exactly
# the assemblage's sites, in the same first-appearance order ComMatrix will use.
sites = unique(phylocom.site)
idx = indexin(sites, string.(env.ID_env))
coord = DataFrame(site = sites,
                  x = (env.xmin[idx] .+ env.xmax[idx]) ./ 2,
                  y = (env.ymin[idx] .+ env.ymax[idx]) ./ 2)

birds = Assemblage(phylocom, coord)

default(color = cgrad(:Spectral, rev = true))
plot(birds)

# birds.nwk stores support values as internal node labels, and support = 1
# recurs thousands of times; Phylo requires unique node names, so strip the
# internal labels (the numeric token right after a close paren) before parsing.
nwk = replace(read(joinpath(datadir, "birds.nwk"), String), r"\)[0-9.]+" => ")")
tree = parsenewick(nwk)
# keep only tips shared with the assemblage (after the Cross.csv relabeling) so
# clade subsetting never references a species that is absent from the data
keeptips!(tree, intersect(getleafnames(tree), speciesnames(birds)))
sort!(tree) # sort the nodes on the tree in order of size - useful for plotting
plot(tree, treetype = :fan, showtips = false, tipfont = (5,))

### Extract information from a single clade

nodes = nodenamefilter(!isleaf, tree)
nodevec = collect(nodes)
first(nodevec, 4)

randnode = nodevec[131]

first(nodespecies(tree, randnode), 4)

rand_clade = get_clade(birds, tree, randnode)
plot(rand_clade, title = randnode)

### Comparing the richness of sister clades

plot_node(birds, tree, randnode)

### Using randomization to assess significance of distribution differences

ch1, ch2 = getchildren(tree, randnode)[1:2]
sims = simulate_descendants(rand_clade, tree, ch1; method = :tipshuffle)
SOS = calculate_SOS(sims)
plot(SOS, rand_clade, clim = (-8,8), fillcolor = :RdYlBu, title = "SOS for clade $randnode")

GND = calculate_GND(sims)

### Putting it all together

# use as
SOS, GND = process_node(birds, tree, randnode; method = :tipshuffle)

### Calculate GND for all nodes

# Scan node divergence across the whole tree. `node_gnd` returns just the per-node
# GND (no SOS maps), and defaults to the fast :tipshuffle null, so it runs over
# the full bird tree in reasonable time. (`node_based_analysis(birds, tree;
# method = :tipshuffle)` is still available if you also want the per-cell SOS maps.)
gnd = node_gnd(birds, tree)   # Dict(nodename => GND), defaults to :tipshuffle

# strongly divergent nodes: GND > 0.8
divergent = [node for (node, g) in gnd if !isnan(g) && g > 0.8]
sort(divergent, by = n -> gnd[n], rev = true)   # inspect: most divergent first

# Map GND onto the tree, showing ONLY the divergent nodes. Phylo's recipe draws a
# marker at every node and (in current Plots) renders NaN-marker_z nodes as solid
# black, so colour can't hide the rest. Instead give each node an explicit size -
# 6 for divergent nodes, 0 (nothing drawn) otherwise - in the recipe's own node
# order, which comes from the layout helper Phylo._findxy.
divset = Set(divergent)
layoutnodes = Phylo._findxy(tree)[3]
markersizes = [node in divset ? 6 : 0 for node in layoutnodes]

plot(tree, treetype = :fan, showtips = false,
     marker_z = Dict(node => gnd[node] for node in divergent),
     markersize = markersizes, markerstrokewidth = 0,
     color = cgrad(:YlOrRd, 10, categorical = true),
     size = (1000, 1000), clim = (0, 1)
     )

### Similarity of SOS patterns among strongly divergent nodes

# Using `divergent` (GND > 0.8, defined above): compute each node's per-cell SOS
# pattern, correlate pairwise (over cells where both are defined - SOS is NaN
# where a clade is absent), convert to distances (1 - |r|), and ordinate with MDS
# so nodes with similar SOS patterns sit close together.

# cells x nodes matrix of SOS patterns
sosmat = reduce(hcat, process_node(birds, tree, node; method = :tipshuffle)[1] for node in divergent)

# pairwise-complete correlations -> distance = 1 - |r|
nnode = length(divergent)
D = zeros(nnode, nnode)
for i in 1:nnode, j in i+1:nnode
    a, b = view(sosmat, :, i), view(sosmat, :, j)
    ok = .!(isnan.(a) .| isnan.(b))
    r = count(ok) > 2 ? cor(a[ok], b[ok]) : 0.0
    D[i, j] = D[j, i] = 1 - abs(isnan(r) ? 0.0 : r)
end

# MDS of the node-distance matrix
nodemds = fit(MDS, D; distances = true, maxoutdim = 2)
nodecoords = predict(nodemds)   # 2 x nnode

plt = scatter(nodecoords[1, :], nodecoords[2, :], label = "",
    xlabel = "MDS axis 1", ylabel = "MDS axis 2",
    title = "Similarity of SOS patterns (nodes with GND > 0.8)")

# put each node label just above its point (offset by a small fraction of the
# y-range, anchored at the text's bottom) so the markers don't cover the labels
yr = extrema(nodecoords[2, :])
dy = 0.025 * (yr[2] - yr[1])
annotate!(plt, [(nodecoords[1, i], nodecoords[2, i] + dy, text(divergent[i], 6, :bottom))
                for i in eachindex(divergent)])


# plot the parent-vs-descendants distributions for the most divergent node.
# (Pick from `divergent` rather than hard-coding a name: Node N labels are assigned
# at parse time and only exist on this exact tree - check with `hasnode(tree, n)`.)
plot_node(birds, tree, "Node 19946")