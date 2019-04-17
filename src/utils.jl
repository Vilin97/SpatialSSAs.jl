using ReactionNetworkImporters, Random, DiffEqProblemLibrary.JumpProblemLibrary, RandomNumbers.Xorshifts, DiffEqJump, Parameters

struct SpatialJumpProb{J, I, F, S, M1, M2}
    jump_prob::J
    box_width::I
    dimension::I
    diff_const::F
    stationary_species::S
    initial_state::M1
    connectivity_matrix::M2
end

function rx_update_props!(vol_rates, vol_rx_rates, vol_diff_rates, rx_matrix, diff_matrix, massaction_jump, dep_graph, connectivity_matrix, j, mu, state, diff_rate)
    j_state = state[j]
    for k in dep_graph[mu]
        old_prop = rx_matrix[j][k]
        prop = DiffEqJump.evalrxrate(j_state, k, massaction_jump)
        rx_matrix[j][k] = prop
        vol_rx_rates[j] += (prop - old_prop)
    end
    L = Float64(length(connectivity_matrix[j]))
    for (n,_) in massaction_jump.net_stoch[mu]
        old_prop = diff_matrix[j][n]
        prop = L * diff_rate * j_state[n]
        diff_matrix[j][n] = prop
        vol_diff_rates[j] += (prop - old_prop)
    end
    vol_rates[j] = vol_rx_rates[j] + vol_diff_rates[j]
end

function diff_update_props!(vol_rates, vol_rx_rates, vol_diff_rates, rx_matrix, diff_matrix, massaction_jump, connectivity_matrix, j, d, n, state, diff_rate, spec_to_dep_rxs,L_vecs)
    for k in spec_to_dep_rxs[n]
        old_prop = rx_matrix[j][k]
        prop = DiffEqJump.evalrxrate(state[j], k, massaction_jump)
        rx_matrix[j][k] = prop
        vol_rx_rates[j] += (prop - old_prop)

        old_prop = rx_matrix[d][k]
        prop = DiffEqJump.evalrxrate(state[d], k, massaction_jump)
        rx_matrix[d][k] = prop
        vol_rx_rates[d] += (prop - old_prop)
    end

    L = L_vecs[length(connectivity_matrix[j])]
    old_prop = diff_matrix[j][n]
    prop = L * diff_rate * state[j][n]
    diff_matrix[j][n] = prop
    vol_diff_rates[j] += (prop - old_prop)
    vol_rates[j] = vol_rx_rates[j] + vol_diff_rates[j]

    L = L_vecs[length(connectivity_matrix[d])]
    old_prop = diff_matrix[d][n]
    prop = L * diff_rate * state[d][n]
    diff_matrix[d][n] = prop
    vol_diff_rates[d] += (prop - old_prop)
    vol_rates[d] = vol_rx_rates[d] + vol_diff_rates[d]
end
function get_rx_props(state, j, massaction_jump, M)
    j_state = state[j]
    return [DiffEqJump.evalrxrate(j_state, k, massaction_jump) for k in 1:M]
end

function get_diff_props(state, j, diff_rate, connectivity_matrix)
    j_state = state[j]
    L = Float64(length(connectivity_matrix[j]))
    diff_props = (L * diff_rate) .* j_state
    return diff_props
end

function record_state!(states, save_all_jumps,state,times,t,idx)
    if save_all_jumps
        push!(states, copy(state))
        push!(times, t)
    else
        temp = idx[1]
        while temp <= length(times) && times[temp] < t
            states[temp] = copy(state)
            temp += 1
        end
        idx[1] = temp
    end
end

function fill_states!(states,save_all_jumps,state,times,idx)
    if save_all_jumps
        temp = idx[1]
        while temp < length(times)
            states[idx] = state
            temp += 1
        end
    end
end

function DM_search(total_rate, propensities, s)
    r = total_rate*rand(s)
    run_sum = 0.0
    for k in 1:length(propensities)
        run_sum += propensities[k]
        if run_sum > r
            return k
        end
    end
    return 0
end

phi(j,m) = ( (j-1)%m+1,(div(j-1,m))%m+1,(div(j-1,m^2)+1) )
phi_inverse(x,y,z,m) = x + (y-1)*m + (z-1)*m^2

function get_connectivity_matrix(box_width, dimension)
    vol_num = box_width^dimension
    connectivity_matrix = Array{Int,1}[]
    for j in 1:vol_num
        x,y,z = phi(j,box_width)
        if dimension == 1
            potential_neighbors = [(x-1,y,z), (x+1,y,z)]
        elseif dimension == 2
            potential_neighbors = [(x-1,y,z), (x+1,y,z), (x,y-1,z), (x,y+1,z)]
        elseif dimension == 3
            potential_neighbors = [(x-1,y,z), (x+1,y,z), (x,y-1,z), (x,y+1,z), (x,y,z-1), (x,y,z+1)]
        end
        real_neighbors = Int[]
        for (x,y,z) in potential_neighbors

            if 1<=x<=box_width && 1<=y<=box_width && 1<=z<=box_width
                push!(real_neighbors, phi_inverse(x,y,z,box_width))
            end
        end
        #println("j = $j. neighbors are: $real_neighbors")
        push!(connectivity_matrix, real_neighbors)
    end
    return connectivity_matrix
end

function get_spatial_jump_prob(jump_prob, box_width, dimension, diff_const, stationary_species, initial_state)
    connectivity_matrix = get_connectivity_matrix(box_width, dimension)
    return SpatialJumpProb(jump_prob, box_width, dimension, diff_const, stationary_species, initial_state, connectivity_matrix)
end

#load specific spatial_jump_probs
function dna_spatial(box_width = 21, dimension = 1, diff_const = 10^-4)
    jump_prob = dna_jump_prob(Direct())
    stationary_species = Set([1,4])
    vol_num = box_width^dimension
    N = length(jump_prob.prob.u0)
    initial_state = [zeros(Int,N) for i in 1:vol_num]
    c = div(vol_num,2)+1
    initial_state[c] = jump_prob.prob.u0
    return get_spatial_jump_prob(jump_prob, box_width, dimension, diff_const, stationary_species, initial_state)
end

# functions to load specific models

function dna_jump_prob(method)
    JumpProblemLibrary.importjumpproblems()
    prob = prob_jump_dnarepressor
    return JumpProblem(prob.discrete_prob,method,prob.network, save_positions = (false,false))
end

function multi_jump_prob(method)
    JumpProblemLibrary.importjumpproblems()
    prob = prob_jump_multistate
    return JumpProblem(prob.discrete_prob,method,prob.network, save_positions = (false,false))
end

function load_model(model_name::String, method)
    if model_name == "DNA"
        jump_prob = dna_jump_prob(method)
    elseif model_name == "B_cell"
        jump_prob = big_jump_prob(method)
    elseif model_name == "multistate"
        jump_prob = multi_jump_prob(method)
    end
    return jump_prob
end

#plotting

function plot_nth_frame(n, k, m, times, sstates)
    max_val_1=maximum([maximum(sstates[n][:,k]) for n in 1:num_of_sample_times])
    max_val_2=maximum([maximum(sstates[n][:,m]) for n in 1:num_of_sample_times])
    p1 = bar(sstates[n][:,k], label=prob.network.syms[k],ylim=(0,max_val_1))
    p2 = bar(sstates[n][:,m], label=prob.network.syms[m],ylim=(0,max_val_2))
    plot(p1,p2,layout=(2,1))
end
