# construct the extensive formulation

using JuMP, Gurobi
# using Statistics, CSV, DataFrames;
# using Random;

function find_topology(k, fs, gs)

    ### find network topology given degree sequence, steady states, and network dynamics
    # xs: list of steady states; k: list of degree sequence; fs, gs: function values for self dynamics and interaction between adjacent nodes at steady states
    ###

    # initialize model
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(), "OutputFlag" => 0));

    # define variables
    n_node = length(k);
    @variable(mp, Adj[i in 1:n_node, j in 1:n_node], Bin);
    @variable(mp, aux_var[i in 1:n_node]);

    # constraints
    # symmetry
    @constraint(mp, symmetry[i in 1:n_node, j in 1:n_node], Adj[i,j] == Adj[j,i]);
    # degree sequence
    @constraint(mp, degree[i in 1:n_node], sum(Adj[i,j] for j in 1:n_node) == k[i]);

    # objective function
    @expression(mp, cost[i in 1:n_node], fs[i] + sum(Adj[i,j]*gs[i,j] for j in 1:n_node) );
    # linearize the absolute term in the obj
    @constraint(mp, LB1[i in 1:n_node], aux_var[i] >= cost[i]);
    @constraint(mp, LB2[i in 1:n_node], aux_var[i] >= -cost[i]);

    @objective(mp, Min, sum(aux_var[i] for i in 1:n_node) );

    t0 = time();
    optimize!(mp);
    t1 = time();

    t_sol = t1 - t0

    Adj = mp[:Adj];

    return Adj, t_sol;
end

# test




