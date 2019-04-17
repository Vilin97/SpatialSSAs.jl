module SpatialSSAs

export my_solve, get_rx_props, get_diff_props

include("utils.jl")
#include("RSSA.jl")
#include("CR_Structures.jl")
#include("CR.jl")
include("Direct_NRM.jl")
#include("NRM.jl")
#include("benchmark.jl")
#include("test.jl")

function my_solve(spatial_jump_prob, method::String, times=nothing)
    if method == "Direct"
        times, states = dm_nrm_solve(spatial_jump_prob, "Direct", times=times)
    elseif method == "NRM"
        times, states = dm_nrm_solve(spatial_jump_prob, "NRM", times=times)
    end
    return times, states
end

end # module
