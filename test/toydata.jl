# A small data set for testing the analysis end to end: twelve species in two clades of
# six, on a 6 x 5 grid. Clade X lives mostly in the western half and clade Y mostly in the
# eastern half, so the root is divergent. Site 1 holds every species, so the null model
# cannot vary there for any node, and site 30 is empty.
#
# The analysable nodes are those whose two children each have at least three species:
# the root, X and Y.
const TOY_NEWICK =
    "((((a:1,b:1):1,c:2):1,((d:1,e:1):1,f:2):1)X:1," *
    "(((g:1,h:1):1,i:2):1,((j:1,k:1):1,l:2):1)Y:1)root;"
const TOY_ANALYSABLE = ["root", "X", "Y"]

function toy_data()
    rng = Xoshiro(1)
    coords = Float64[repeat(1:6, 5) repeat(1:5; inner=6)]
    west = coords[:, 1] .<= 3
    occ = falses(12, 30)
    for species in 1:12, site in 1:30
        home = (species <= 6) == west[site]
        occ[species, site] = rand(rng) < (home ? 0.7 : 0.15)
    end
    occ[:, 1] .= true
    occ[:, 30] .= false
    sites = ["s$i" for i in 1:30]
    assemblage = Assemblage(Matrix(occ), coords, sites, string.('a':'l'))
    return assemblage, parsenewick(TOY_NEWICK)
end
