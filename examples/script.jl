using CSV, DataFrames, SpatialEcology, Phylo, Plots
using Distances, MultivariateStats

using Nodiv

# Data files live in the local SpatialEcology dev copy's docs.
const datadir = "/Users/cvg147/.julia/dev/SpatialEcology/docs/data"

### Load data and create objects

phylocom = CSV.read(joinpath(datadir, "tyrann_phylocom.tsv"), DataFrame)
first(phylocom, 4)

coord = CSV.read(joinpath(datadir, "tyrann_coords.tsv"), DataFrame)
first(coord, 4)

phylocom.Plot = string.(phylocom.Plot)
coord.cell = string.(coord.cell)
tyrants = Assemblage(phylocom, coord)

default(color = cgrad(:Spectral, rev = true))
plot(tyrants)

### Ordinate sites in species space
#
# NOTE: MultivariateStats provides *metric* MDS (classical / PCoA-style), not
# non-metric NMDS — there is no off-the-shelf Kruskal NMDS in the Julia
# ecosystem. This is the closest readily available ordination; swap in a
# dedicated NMDS implementation if strict non-metric scaling is required.

# Jaccard dissimilarity between sites (presence/absence data). `occurrences`
# returns a species-by-site matrix, so sites are the columns (dims = 2).
# Qualify `pairwise` because both Distances and SpatialEcology export it.
sitedist = Distances.pairwise(Jaccard(), Matrix(occurrences(tyrants)); dims = 2)

mds = fit(MDS, sitedist; distances = true, maxoutdim = 2)
mdscoords = predict(mds)  # 2 x nsites

scatter(mdscoords[1, :], mdscoords[2, :],
    marker_z = richness(tyrants), label = "",
    xlabel = "MDS axis 1", ylabel = "MDS axis 2",
    colorbar_title = "richness", title = "Site ordination (metric MDS, Jaccard)")

tree = open(parsenewick, joinpath(datadir, "tyrannid_tree.tre"))
sort!(tree) # sort the nodes on the tree in order of size - useful for plotting
plot(tree, treetype = :fan, tipfont = (5,))

### Extract information from a single clade

nodes = nodenamefilter(!isleaf, tree)
nodevec = collect(nodes)
first(nodevec, 4)

randnode = nodevec[131]

first(nodespecies(tree, randnode), 4)

rand_clade = get_clade(tyrants, tree, randnode)
plot(rand_clade, title = randnode)

### Comparing the richness of sister clades

plot_node(tyrants, tree, randnode)

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
SOS, GND = process_node(tyrants, tree, randnode)

### Final step: Applying the method to the entire tree

SOSs, GNDs = node_based_analysis(tyrants, tree);

plot(tree,
     showtips = false, marker_z = GNDs,
     color = cgrad(:YlOrRd, 10, categorical = true),
     markersize = 15 .* GNDs, markerstrokewidth = 0,
     size = (600, 1000), clim = (0,1)
     )
