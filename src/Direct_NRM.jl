using Random, Parameters, RandomNumbers.Xorshifts,  DiffEqBiological, DiffEqJump, DataStructures

function dm_nrm_solve(spatial_jump_prob, method="Direct"; times=nothing)
    parameters = dm_nrm_do_setup(spatial_jump_prob, method, times)
    return dm_nrm_run_simulation(parameters...)
end

function dm_nrm_do_setup(spatial_jump_prob, method, times)
    @unpack jump_prob, initial_state = spatial_jump_prob
    @unpack massaction_jump = jump_prob
    @unpack u0 = jump_prob.prob
    t0 = jump_prob.prob.tspan[1]
    N = length(u0)
    dep_graph = DiffEqJump.make_dependency_graph(N, massaction_jump)
    if times === nothing
        times = [t0]
        states = [initial_state]
        save_all_jumps = true
    else
        states = [initial_state for i in times]
        save_all_jumps = false
    end
    return spatial_jump_prob, method, dep_graph, times, states, save_all_jumps
end


function dm_nrm_run_simulation(spatial_jump_prob, method, dep_graph, times, states, save_all_jumps)
    @unpack jump_prob, box_width, dimension, diff_const, stationary_species, initial_state, connectivity_matrix =
    spatial_jump_prob
    diff_rate = diff_const*box_width^2
    vol_num = box_width^dimension
    t0,t_end = jump_prob.prob.tspan
    @unpack massaction_jump = jump_prob
    M = get_num_majumps(massaction_jump)
    N = length(jump_prob.prob.u0)
    spec_to_dep_rxs = DiffEqJump.spec_to_dep_rxs_map(N, massaction_jump)
    state = deepcopy(initial_state)
    rx_matrix = [zeros(Float64, M) for i in 1:vol_num]
    diff_matrix = [zeros(Float64, N) for i in 1:vol_num]
    vol_rx_rates = zeros(vol_num)
    vol_diff_rates = zeros(vol_num)
    vol_rates = zeros(vol_num)
    for j in 1:vol_num
        rx_matrix[j] = get_rx_props(state, j, massaction_jump, M)
        diff_matrix[j] = get_diff_props(state, j, diff_rate, connectivity_matrix)
        vol_rx_rates[j] = sum(rx_matrix[j])
        vol_diff_rates[j] = sum(diff_matrix[j])
        vol_rates[j] = vol_rx_rates[j] + vol_diff_rates[j]
    end
    total_rate = sum(vol_rates)
    if method == "NRM"
        h = MutableBinaryMinHeap(Random.randexp(vol_num)./vol_rates) end
    t = t0
    rx_counter = 0
    diff_counter = 0
    idx = [2]
    ones_vecs = [ones(k) for k in 1:2*dimension]
    L_vecs = [Float64(L) for L in 1:2*dimension]
    s = Xoroshiro128Plus()
    while t < t_end
        if method == "Direct"
            t += Random.randexp(s)/total_rate
            j = DM_search(total_rate, vol_rates, s)
        elseif method == "NRM"
            t, j = top_with_handle(h)
        end
        rx_rate = vol_rx_rates[j]
        rx_occurs = rand(s)*vol_rates[j] < rx_rate #rx occur
        if rx_occurs
            rx_counter += 1
            propensities = rx_matrix[j]
            mu = DM_search(rx_rate, propensities, s)
            j_state = state[j]
            DiffEqJump.executerx!(j_state, mu, massaction_jump)
            state[j] = j_state
            total_rate -= vol_rates[j]
            rx_update_props!(vol_rates, vol_rx_rates, vol_diff_rates, rx_matrix, diff_matrix, massaction_jump, dep_graph, connectivity_matrix, j, mu, state, diff_rate)
            total_rate += vol_rates[j]
            if method == "NRM"
                if vol_rates[j] > 0 update!(h, j, t+Random.randexp(s)/vol_rates[j])
                else update!(h, j, Inf) end
            end
        else #diffusion occurs
            diff_prop = vol_diff_rates[j]
            n = DM_search(diff_prop, diff_matrix[j], s)
            if n in stationary_species
                continue
            end
            diff_counter += 1
            neighbors = connectivity_matrix[j]
            L = length(neighbors)
            temp = DM_search(L_vecs[L], ones_vecs[L], s)
            d = neighbors[temp]
            state[j][n] -= 1
            state[d][n] += 1
            total_rate -= vol_rates[j]
            total_rate -= vol_rates[d]
            diff_update_props!(vol_rates, vol_rx_rates, vol_diff_rates, rx_matrix, diff_matrix, massaction_jump, connectivity_matrix, j, d, n, state, diff_rate, spec_to_dep_rxs,L_vecs)
            total_rate += vol_rates[j]
            total_rate += vol_rates[d]
            if method == "NRM"
                if vol_rates[j] > 0 update!(h, j, t+Random.randexp(s)/vol_rates[j])
                else update!(h, j, Inf) end
                if vol_rates[d] > 0 update!(h, d, t+Random.randexp(s)/vol_rates[d])
                else update!(h, d, Inf) end
            end
        end
        record_state!(states, save_all_jumps, state, times, t, idx)
    end
    println("number of reactions: $rx_counter. number of diffusions: $diff_counter")
    fill_states!(states,save_all_jumps,state,times,idx)
    return times, states
end
