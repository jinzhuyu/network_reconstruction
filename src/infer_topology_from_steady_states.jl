# construct the extensive formulation

using JuMP, Gurobi
# using Statistics, CSV, DataFrames;
# using Random;

function find_topology(k, fs, gs)

    ### find network topology given degree sequence, steady states, and network dynamics
    # ks: list of degree sequence; fs, gs: function values for self dynamics and interaction between adjacent nodes at steady states
    ###

    # initialize model
    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(), "OutputFlag" => 0));

    # define variables
    n_node = length(k);
    @variable(model, Adj[i in 1:n_node, j in 1:n_node], Bin);
    @variable(model, aux_var[i in 1:n_node]);

    # constraints
    # symmetry
    @constraint(model, symmetry[i in 1:n_node, j in 1:n_node], Adj[i,j] == Adj[j,i]);
    # degree sequence
    @constraint(model, degree[i in 1:n_node], sum(Adj[i,j] for j in 1:n_node) == k[i]);

    # objective function
    @expression(model, cost[i in 1:n_node], fs[i] + sum(Adj[i,j]*gs[i,j] for j in 1:n_node) );
    # linearize the absolute term in the obj
    @constraint(model, LB1[i in 1:n_node], aux_var[i] >= cost[i]);
    @constraint(model, LB2[i in 1:n_node], aux_var[i] >= -cost[i]);

    @objective(model, Min, sum(aux_var[i] for i in 1:n_node) );

    # solve the model
    t0 = time();
    optimize!(model);
    t1 = time();

    if termination_status(m) != MOI.OPTIMAL
        @warn("===== Final solution not optimal =====")
        return
    end

    # return solution time and solution
    t_sol = t1 - t0
    Adj = model[:Adj];

    return Adj, t_sol;
end

# test




