module Nodiv

include("functions.jl")
# plot_node / plot_node! are exported by the @userplot macro in functions.jl
export nodespecies, get_clade
export simulate_descendants
export calculate_GND, calculate_SOS
export process_node, node_based_analysis, node_gnd

end
