using CSV, DataFrames, SpatialEcology, Phylo, Plots
using Distances, MultivariateStats

using Nodiv

const datadir = "/Users/cvg147/Dropbox/Arbejde/Current projects/NodivWorkshop"

### Load data and create objects

pam = CSV.read(joinpath(datadir, "PAM_E.csv"), DataFrame)
first(pam, 4)

env = CSV.read(joinpath(datadir, "Env.csv"), DataFrame)
first(env, 4)

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

### Ordinate sites in species space
#
# NOTE: MultivariateStats provides *metric* MDS (classical / PCoA-style), not
# non-metric NMDS — there is no off-the-shelf Kruskal NMDS in the Julia
# ecosystem. This is the closest readily available ordination; swap in a
# dedicated NMDS implementation if strict non-metric scaling is required.

# Jaccard dissimilarity between sites (presence/absence data). `occurrences`
# returns a species-by-site matrix, so sites are the columns (dims = 2).
# Qualify `pairwise` because both Distances and SpatialEcology export it.
sitedist = Distances.pairwise(Jaccard(), Matrix(occurrences(birds)); dims = 2)

mds = fit(MDS, sitedist; distances = true, maxoutdim = 2)
mdscoords = predict(mds)  # 2 x nsites

scatter(mdscoords[1, :], mdscoords[2, :],
    marker_z = richness(birds), label = "",
    xlabel = "MDS axis 1", ylabel = "MDS axis 2",
    colorbar_title = "richness", title = "Site ordination (metric MDS, Jaccard)")

# birds.nwk stores support values as internal node labels, and support = 1
# recurs thousands of times; Phylo requires unique node names, so strip the
# internal labels (the numeric token right after a close paren) before parsing.
nwk = replace(read(joinpath(datadir, "birds.nwk"), String), r"\)[0-9.]+" => ")")
tree = parsenewick(nwk)
# keep only tips shared with the assemblage (after the Cross.csv relabeling) so
# clade subsetting never references a species that is absent from the data
keeptips!(tree, intersect(getleafnames(tree), speciesnames(birds)))
sort!(tree) # sort the nodes on the tree in order of size - useful for plotting
plot(tree, treetype = :fan, tipfont = (5,))

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

rmg = matrixrandomizer(rand_clade)
newcomm = rand(rmg)
plot(newcomm, title = "randomized version of $randnode")

ch1, ch2 = getchildren(tree, randnode)[1:2]
plot_node(newcomm, tree, randnode)

sims = simulate_descendants(rand_clade, tree, ch1)
SOS = calculate_SOS(sims)
plot(SOS, rand_clade, clim = (-8,8), fillcolor = :RdYlBu, title = "SOS for clade $randnode")

GND = calculate_GND(sims)
first(GND, 4)

### Putting it all together

# use as
SOS, GND = process_node(birds, tree, randnode)

### Calculate GND for all nodes

# Scan node divergence across the whole tree. `node_gnd` returns just the per-node
# GND (no SOS maps), and defaults to the fast :tipshuffle null, so it runs over
# the full bird tree in reasonable time. (`node_based_analysis(birds, tree;
# method = :tipshuffle)` is still available if you also want the per-cell SOS maps.)
nodes, GNDs = node_gnd(birds, tree)

# the most divergent nodes
valid = findall(!isnan, GNDs)
ranked = valid[sortperm(GNDs[valid], rev = true)]
[nodes[ranked] GNDs[ranked]][1:10, :]

# map GND onto the tree
plot(tree,
     showtips = false, marker_z = GNDs,
     color = cgrad(:YlOrRd, 10, categorical = true),
     markersize = 15 .* GNDs, markerstrokewidth = 0,
     size = (600, 1000), clim = (0,1)
     )
