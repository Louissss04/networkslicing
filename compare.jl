
using JuMP
using CPLEX
using Gurobi
using Printf
using DelimitedFiles
using Statistics
import MathOptInterface
const MOI = MathOptInterface
using Random

include("ReadData.jl")
include("constraint.jl")

Profun = 3
nP = 2
ntests = 3

graphFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-1-7-graph_info-" * string(ntests)
flowFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-1-7-flow_info-" * string(ntests)

# readFlowFile(flowFile)

# 调用 readProb 函数并传递 graphFile、flowFile 和 ntests 参数
nI, nL, nV, nF, L, V, Link_cap, Link_reliability, Link_delay, Node_cap, Node_delay, Node_reliability, outgoing_cap, incoming_cap, V_processfunction, nK, nS, PairC, DataRate, FunctionChainC, E2EDelay, E2Ereliability = readProb(graphFile, flowFile, ntests)

# create model
# model = Model(CPLEX.Optimizer)
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "MIPGap", 0.005)


# def variable
@variable(model, x[1:nV, 1:nK, 1:nS], Bin)  # x[v, k, s]
@variable(model, x0[1:nV, 1:nK], Bin)  # x[v, k]
@variable(model, y[1:nV], Bin)  # y[v]
@variable(model, theta[1:nK, 1:nS+1] >= 0) # theta[k, s]

@variable(model, theta_L[1:nK] >= 0) # theta[L, k]
@variable(model, theta_N[1:nK] >= 0) # theta[N, k]
@variable(model, rel_L[1:nK])   #
@variable(model, rel_N[1:nK])   #

@variable(model, lambda[1, 1:nK, 1:nS+1, 1:nP] >= 0) # r[k, s, p]
@variable(model, z[1:nL, 1:nK, 1:nS+1, 1:nP], Bin)   # z[L, k, s, p]
@variable(model, z0[1:nL, 1:nK], Bin)                # z[L, k]
@variable(model, r[1:nL, 1:nK, 1:nS+1, 1:nP] >= 0)   # r[L, k, s, p]
@variable(model, R[1:nL, 1:nK] >= 0)                 # R[L, k]
@variable(model, R0[1:nV, 1:nK] >= 0)                # R[v, k]

# create constraint
VNFplacement(model, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0)

#TrafficroutingMINLP(model, x, r, z, lambda, V, L, PairC, nV, nI, nL, nK, nS, nP, z0)

TrafficroutingMILP(model, x, r, z, V, L, PairC, nV, nI, nL, nK, nS, nP, z0)

Linkcap(model, r, R, DataRate, Link_cap, nL, nK)

E2EreliabilityCons(model, x, x0, r, z, z0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, nP)

E2EdelayCons(model, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP)

# add bound constraint
@constraint(model, r .<= 1)
@constraint(model, lambda .<= 1)


# objective function

#@objective(model, Min, sum(y))

@objective(model, Min, sum(y)+0.005sum(r))

# solve
optimize!(model)

objective_value1 = JuMP.objective_value(model)

# create model
# model = Model(CPLEX.Optimizer)
model0 = Model(Gurobi.Optimizer)

# def variable
@variable(model0, x[1:nV, 1:nK, 1:nS], Bin)  # x[v, k, s]
@variable(model0, x0[1:nV, 1:nK], Bin)  # x[v, k]
@variable(model0, y[1:nV], Bin)  # y[v]
@variable(model0, theta[1:nK, 1:nS+1] >= 0) # theta

@variable(model0, theta_L[1:nK] >= 0) # theta[L, k]
@variable(model0, theta_N[1:nK] >= 0) # theta[N, k]
@variable(model0, rel_L[1:nK])   #
@variable(model0, rel_N[1:nK])   #

@variable(model0, lambda[1, 1:nK, 1:nS+1, 1:nP] >= 0) # r[k, s, p]
@variable(model0, z[1:nL, 1:nK, 1:nS+1, 1:nP], Bin)   # z[L, k, s, p]
@variable(model0, z0[1:nL, 1:nK], Bin)                # z[L, k]
@variable(model0, r[1:nL, 1:nK, 1:nS+1, 1:nP] >= 0)   # r[L, k, s, p]
@variable(model0, R[1:nL, 1:nK] >= 0)                 # R[L, k]
@variable(model0, R0[1:nV, 1:nK] >= 0)                # R[v, k]

# create constraint
VNFplacement(model0, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0)

TrafficroutingMINLP(model0, x, r, z, lambda, V, L, PairC, nV, nI, nL, nK, nS, nP, z0)

#TrafficroutingMILP(model0, x, r, z, V, L, PairC, nV, nI, nL, nK, nS, nP, z0)

Linkcap(model0, r, R, DataRate, Link_cap, nL, nK)

E2EreliabilityCons(model0, x, x0, r, z, z0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, nP)

E2EdelayConsforMINLP(model0, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP)

# add bound constraint
@constraint(model0, r .<= 1)
@constraint(model0, lambda .<= 1)


# objective function

#@objective(model, Min, sum(y))

@objective(model0, Min, sum(y)+0.005sum(r))

# solve
optimize!(model0)

objective_value2 = JuMP.objective_value(model0)

println("MINLP Objective Value: $objective_value1")
println("MILP Objective Value: $objective_value2")

compare(model, model0)