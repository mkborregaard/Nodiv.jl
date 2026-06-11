# Node-based analysis of bird diversity, run in parallel in two spaces:
#   birds_e - environmental (PC1/PC2) space  (PAM_E + Env.csv)
#   birds_g - geographic space               (PAM_G + g_space shapefile)
# Reusable, general steps live in the Nodiv package; everything dataset-specific
# (file layout, taxonomy crosswalk, coordinate building) stays in this script.

using CSV, DataFrames, SpatialEcology, Phylo, Plots, Shapefile
using MultivariateStats, Statistics, JLD2
using Nodiv

datadir = "/Users/cvg147/Dropbox/Arbejde/Current projects/NodivWorkshop"
default(color = cgrad(:Spectral, rev = true))

### Dataset-specific helpers ---------------------------------------------------

# Cross.csv maps the PAM/BirdLife taxonomy (Species1) to the tree/BirdTree
# taxonomy (Species3); relabel PAM species to the tree names so the datasets
# share as many taxa as possible.
cross = CSV.read(joinpath(datadir, "Cross.csv"), DataFrame)
namemap = Dict(string(r.Species1) => replace(string(r.Species3), " " => "_")
               for r in eachrow(cross) if !ismissing(r.Species1) && !ismissing(r.Species3))

# read a Newick tree, stripping numeric internal (support) labels that Phylo
# rejects as duplicate node names
readtree(path) = parsenewick(replace(read(path, String), r"\)[0-9.]+" => ")"))

# build an Assemblage from long-format occurrences + a (site, x, y) coordinate
# lookup: relabel species to the tree taxonomy, and align the coordinates to the
# matrix's sites (SpatialEcology matches by row order, not by name).
function build_assemblage(species_raw, sitevals, coords)
    species = [get(namemap, s, replace(s, " " => "_")) for s in species_raw]
    pc = DataFrame(site = string.(sitevals), abundance = 1, species = species)
    sites = unique(pc.site)
    idx = indexin(sites, string.(coords.site))
    Assemblage(pc, DataFrame(site = sites, x = coords.x[idx], y = coords.y[idx]))
end

### Environmental-space assemblage ---------------------------------------------

env = CSV.read(joinpath(datadir, "Env.csv"), DataFrame)
pam_e = CSV.read(joinpath(datadir, "PAM_E.csv"), DataFrame)
# env-space coordinates: each cell's PC bin midpoint
coords_e = DataFrame(site = string.(env.ID_env),
                     x = (env.xmin .+ env.xmax) ./ 2,
                     y = (env.ymin .+ env.ymax) ./ 2)
birds_e = build_assemblage(pam_e.Species, pam_e.ID_env, coords_e)
addsitestats!(birds_e, env, :ID_env)        # attach PC bins, area, occupancy, ...
plot(birds_e)

### Geographic-space assemblage ------------------------------------------------

shp = Shapefile.Table(joinpath(datadir, "g_space", "BehrmannMeterGrid_WGS84_land_PCA_30.shp"))
centroid(g) = (ex = extrema(p.x for p in g.points); ey = extrema(p.y for p in g.points);
               ((ex[1] + ex[2]) / 2, (ey[1] + ey[2]) / 2))
cents = centroid.(Shapefile.shapes(shp))
# The equal-area cells become irregular once reprojected to lon/lat, so snap them
# to a regular grid: bin longitude into ~1-degree columns and sin(latitude) into
# rows (which makes the equal-area rows evenly spaced), then map each bin to a
# contiguous integer index so SpatialEcology builds a clean GridData (gridvar
# needs uniform spacing). Warps the map but keeps every cell; a few cells share a
# grid square - that only matters for the image, not the analysis.
# (returns Float indices - SpatialEcology's grid indexing requires float coords)
gridindex(v) = (u = sort(unique(v)); pos = Dict(u .=> eachindex(u)); Float64[pos[x] for x in v])
coords_g = DataFrame(site = string.(shp.ID_geo),
                     x = gridindex(round.(Int, first.(cents))),
                     y = gridindex(round.(Int, sin.(deg2rad.(last.(cents))) .* 180)))
geo_attrs = select(DataFrame(shp), Not(:geometry))
geo_attrs.ID_geo = string.(geo_attrs.ID_geo)

pam_g = CSV.read(joinpath(datadir, "PAM_G.csv"), DataFrame)
birds_g = build_assemblage(pam_g.Species, pam_g.ID_geo, coords_g)
addsitestats!(birds_g, geo_attrs, :ID_geo)  # attach CHELSA bioclim, PC1-3, area, ...
plot(birds_g)

### Shared phylogeny -----------------------------------------------------------

tree = readtree(joinpath(datadir, "birds.nwk"))
prune_to_shared!(tree, birds_e, birds_g)    # keep only taxa present in both
sort!(tree)

### Heavy step: GND + SOS for every node, both spaces, cached to disk ----------
# `node_analysis` computes GND and the per-cell SOS together (SOS is needed for
# GND anyway). This is the slow part - randomisations over the whole tree, and
# ~18k cells for the geographic scan - so cache it: re-running the script just
# reloads the results and jumps straight to the plotting below.
cachefile = joinpath(@__DIR__, "node_analysis.jld2")
if !isfile(cachefile)
    res_e = node_analysis(birds_e, tree)
    res_g = node_analysis(birds_g, tree)
    jldsave(cachefile; res_e, res_g)
end
res_e, res_g = load(cachefile, "res_e", "res_g")   # each a NodeAnalysis (gnd + sos)

### ---- Exploratory plotting (from the cached NodeAnalysis; `_e` vs `_g`) ----- ###

# GND mapped onto the tree (pass the NodeAnalysis, or a filtered Dict for a subset)
plot_gnd(tree, res_e)
plot_gnd(tree, res_g)

# strongly divergent nodes in each space
divergent_e = divergent_nodes(res_e; threshold = 0.8)
divergent_g = divergent_nodes(res_g; threshold = 0.8)

# SOS of the most divergent node mapped onto each space (cached SOS, no recompute)
focal_e = argmax(n -> res_e.gnd[n], divergent_e)
plot(res_e.sos[focal_e], birds_e, fillcolor = :RdYlBu, clim = (-8, 8), title = "env SOS - $focal_e")
focal_g = argmax(n -> res_g.gnd[n], divergent_g)
plot(res_g.sos[focal_g], birds_g, fillcolor = :RdYlBu, clim = (-8, 8), title = "geo SOS - $focal_g")

# parent/SOS/children panel for that node (4th arg = cached SOS, no recompute)
plot_node(birds_e, tree, focal_e, res_e)
plot_node(birds_g, tree, focal_g, res_g)

# ordinate the divergent nodes by SOS-pattern similarity (cached SOS -> distances
# from Nodiv -> MDS; presentation stays here)
function sos_mds_plot(res, nodes, title)
    coords = predict(fit(MDS, sos_distances(res, nodes); distances = true, maxoutdim = 2))
    scatter(coords[1, :], coords[2, :], label = "",
            series_annotations = text.(nodes, 6, :bottom),
            xlabel = "MDS axis 1", ylabel = "MDS axis 2", title = title)
end
sos_mds_plot(res_e, divergent_e, "SOS-pattern similarity (environmental)")
sos_mds_plot(res_g, divergent_g, "SOS-pattern similarity (geographic)")
