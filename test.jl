
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
#ntests = 1
function test(ntests, K, capacityparam)

    println("?????????$ntests???????????")
    file_path = "/share/home/sunqi/GitPublish/ns-algorithm/julia_code/log/MILP-$K-$capacityparam-$ntests.txt"
    file = open(file_path, "a")
    graphFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-$K-$capacityparam-graph_info-$ntests"
    flowFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-$K-$capacityparam-flow_info-$ntests"

    # readFlowFile(flowFile)

    # 调用 readProb 函数并传递 graphFile、flowFile 和 ntests 参数
    nI, nL, nV, nF, L, V, Link_cap, Link_reliability, Link_delay, Node_cap, Node_delay, Node_reliability, outgoing_cap, incoming_cap, V_processfunction, nK, nS, PairC, DataRate, FunctionChainC, E2EDelay, E2Ereliability = readProb(graphFile, flowFile, ntests)

    #println("nK=$nK")
    # create model
    # model = Model(CPLEX.Optimizer)
    model = Model(Gurobi.Optimizer)

    set_optimizer_attribute(model, "Threads", 1)
    set_optimizer_attribute(model, "MIPGap", 0)
    # def variable
    @variable(model, x[1:nV, 1:nK, 1:nS], Bin)  # x[v, k, s]
    @variable(model, x0[1:nV, 1:nK], Bin)  # x[v, k]
    @variable(model, y[1:nV], Bin)  # y[v]
    @variable(model, theta[1:nK, 1:nS+1] >= 0) # theta

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

    #E2EdelayConsforMINLP(model, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP)


    # add bound constraint

    @constraint(model, r .<= 1)
    @constraint(model, lambda .<= 1)


    # objective function

    #@objective(model, Min, sum(y))

    @objective(model, Min, sum(y)+0.005sum(R))

    # solve
    optimize!(model)

    status = JuMP.termination_status(model)

    if status == MOI.OPTIMAL
        rel_total = 0
        for k = 1:nK
            rel = 1
            for v = 1:nV
                if sum(value.(x[v,k,:])) >= 1e-6
                    rel = rel * Node_reliability[v]
                end
            end
            
            for l = 1:nL
                if sum(sum(value.(r[l,k,:,:]))) >= 1e-6
                    rel = rel * Link_reliability[l]
                end
            end
            rel_total += rel
        end

        delay_l = zeros(Float64, nK)
        for k = 1:nK
            for s = 1 : nS+1
                for l = 1:nL
                    if sum(value.(r[l,k,s,:])) >= 1e-6
                        delay_l[k] += Link_delay[l]
                    end
                end
            end
        end


        rel_aver = rel_total / nK
        del_aver = sum(delay_l .+ value.(theta_N)) / nK
        #valuer = value.(r)
        #valuez = value.(z)
        #minus = valuez - valuer
        actived_node = sum(value.(y))
        link_consumption = sum(value.(R))
        reliability = sum(exp.(value.(rel_L) .+ value.(rel_N))) / nK
        delay = sum(value.(theta_L) .+ value.(theta_N)) / nK

        solve_time = MOI.get(model, MOI.SolveTimeSec())

        objective_value = JuMP.objective_value(model)

        #println(file, "minus = $minus")

        println(file, "Objective Value of instance MILP-$K-$capacityparam-$ntests: $objective_value")
        println(file, "Solve time of instance MILP-$K-$capacityparam-$ntests: $solve_time")
        println(file, "Number of actived node of instance MILP-$K-$capacityparam-$ntests: $actived_node")
        println(file, "Total link capacity consumption of instance MILP-$K-$capacityparam-$ntests: $link_consumption")
        println(file, "Average E2E reliability of the services of instance MILP-$K-$capacityparam-$ntests: $reliability")
        println(file, "Real E2E reliability of the services of instance MILP-$K-$capacityparam-$ntests: $rel_aver")
        println(file, "Average E2E delay of the services of instance MILP-$K-$capacityparam-$ntests: $delay")
        println(file, "Real E2E delay of the services of instance MILP-$K-$capacityparam-$ntests: $del_aver")

    else

        print(file, "instance $ntests has No feasible solution!")
        
    end
end

function testlp(ntests, K, capacityparam)

    println("?????????$ntests???????????")
    # file_path = "/share/home/sunqi/GitPublish/ns-algorithm/julia_code/log/MILP-$K-$capacityparam-$ntests.txt"
    # file = open(file_path, "a")
    graphFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-$K-$capacityparam-graph_info-$ntests"
    flowFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-$K-$capacityparam-flow_info-$ntests"

    # readFlowFile(flowFile)

    # 调用 readProb 函数并传递 graphFile、flowFile 和 ntests 参数
    nI, nL, nV, nF, L, V, Link_cap, Link_reliability, Link_delay, Node_cap, Node_delay, Node_reliability, outgoing_cap, incoming_cap, V_processfunction, nK, nS, PairC, DataRate, FunctionChainC, E2EDelay, E2Ereliability = readProb(graphFile, flowFile, ntests)

    #println("nK=$nK")
    # create model
    # model = Model(CPLEX.Optimizer)
    model = Model(Gurobi.Optimizer)

    set_optimizer_attribute(model, "Threads", 1)

    # def variable
    @variable(model, x[1:nV, 1:nK, 1:nS])  # x[v, k, s]
    @variable(model, x0[1:nV, 1:nK])  # x[v, k]
    @variable(model, y[1:nV])  # y[v]
    @variable(model, theta[1:nK, 1:nS+1] >= 0) # theta

    @variable(model, theta_L[1:nK] >= 0) # theta[L, k]
    @variable(model, theta_N[1:nK] >= 0) # theta[N, k]
    @variable(model, rel_L[1:nK])   #
    @variable(model, rel_N[1:nK])   #

    # @variable(model, lambda[1, 1:nK, 1:nS+1, 1:nP] >= 0) # r[k, s, p]
    @variable(model, z[1:nL, 1:nK, 1:nS+1, 1:nP])   # z[L, k, s, p]
    @variable(model, z0[1:nL, 1:nK])                # z[L, k]
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

    #E2EdelayConsforMINLP(model, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP)


    # add bound constraint
    @constraint(model, x .<= 1)
    @constraint(model, x0 .<= 1)
    @constraint(model, y .<= 1)
    @constraint(model, z .<= 1)
    @constraint(model, z0 .<= 1)
    @constraint(model, r .<= 1)
    # @constraint(model, lambda .<= 1)


    # objective function

    #@objective(model, Min, sum(y))

    @objective(model, Min, sum(y)+0.005sum(R))

    # solve
    optimize!(model)

    status = JuMP.termination_status(model)

    if status == MOI.OPTIMAL

        solve_time = MOI.get(model, MOI.SolveTimeSec())

        objective_value = JuMP.objective_value(model)

        println("Objective Value of instance MILP-$K-$capacityparam-$ntests: $objective_value")
        println("Solve time of instance MILP-$K-$capacityparam-$ntests: $solve_time")

    else

        print("instance $ntests has No feasible solution!")
    end
end

# 从命令行获取参数
args = parse.(Int, ARGS)

# 调用 main 函数

test(args...)
