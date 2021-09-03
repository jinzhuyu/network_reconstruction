# construct the extensive formulation

using JuMP, Gurobi;
using NPZ;
# using NLopt;


function find_topology(ks, fs, gs)

    ### find network topology given degree sequence, steady states, and network dynamics
    # ks: list of degree sequence; fs, gs: lists of function values for self dynamics and interaction between adjacent nodes at steady states
    ###

    # initialize model
    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(), "OutputFlag" => 1));

    # define variables
    n_node = length(ks);
    @variable(model, Adj[i in 1:n_node, j in 1:n_node], Bin);
    @variable(model, aux_var[i in 1:n_node]);

    # constraints
    # symmetry
    @constraint(model, symmetry[i in 1:n_node, j in 1:n_node], Adj[i,j] == Adj[j,i]);
    # degree sequence
    @constraint(model, degree[i in 1:n_node], sum(Adj[i,j] for j in 1:n_node) == ks[i]);

    # objective function
    @expression(model, cost[i in 1:n_node], fs[i] + sum(Adj[i,j]*gs[i,j] for j in 1:n_node) );
    # linearize the absolute term in the obj
    @constraint(model, LB1[i in 1:n_node], aux_var[i] >= cost[i]);
    @constraint(model, LB2[i in 1:n_node], aux_var[i] >= -cost[i]);

    @objective(model, Min, sum(aux_var[i] for i in 1:n_node) );

    # @NLobjective(model, Min, sum(cost[i]*cost[i] for i in 1:n_node) );

    # solve the model
    t0 = time();
    optimize!(model);
    t1 = time();

    if termination_status(m) != MOI.OPTIMAL
        @warn("===== Final solution not optimal =====")
        return
    end

    # return solution time and solution
    t_sol = t1 - t0;
    println("===== Solution time for ", n, "-node network is", "t_sol");
    
    obj = objective_value(model);
    println("===== Optimal objective is ", obj);

    Adj = value.(model[:Adj]);

    return obj, Adj, t_sol;
end

# test
function main()
    # load ks, gs, fs
    cd("C:/code/network_reconstruction/data");
    gs = npzread("gs.npy");
    fs = npzread("fs.npy");
    ks = npzread("ks.npy");
    obj, Adj, t_sol = find_topology(ks, fs, gs);
end

main()





