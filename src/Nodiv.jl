module Nodiv

include("functions.jl")
# plot_node / plot_node! are exported by the @userplot macro in functions.jl
export nodespecies, get_clade
export simulate_descendants
export calculate_GND, calculate_SOS
export calculate_GND_rms, calculate_GND_spatial, calculate_GND_ses, gnd_rms, gnd_spatial
export process_node, node_based_analysis, node_analysis, NodeAnalysis
export node_metrics, NodeMetrics
export prune_to_shared!, divergent_nodes, sos_distances

end
