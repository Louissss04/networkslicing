# using JuMP
# using CPLEX
# using Gurobi
# using Printf
# using DelimitedFiles
# using Statistics
# import MathOptInterface
# const MOI = MathOptInterface
# using Random

# include("ReadData.jl")
# include("constraint.jl")

function are_arrays_integers(arrays...; tolerance=1e-5)
    for array in arrays
        for value in array
            if abs(round(value) - value) >= tolerance
                return false
            end
        end
    end
    return true
end
    
    # Profun = 3
    # nP = 2
    
    # #ntests = 1

    # graphFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-40-6-graph_info-" * string(ntests)
    # flowFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-40-6-flow_info-" * string(ntests)

    # # readFlowFile(flowFile)

    # # 调用 readProb 函数并传递 graphFile、flowFile 和 ntests 参数
    # nI, nL, nV, nF, L, V, Link_cap, Link_reliability, Link_delay, Node_cap, Node_delay, Node_reliability, outgoing_cap, incoming_cap, V_processfunction, nK, nS, PairC, DataRate, FunctionChainC, E2EDelay, E2Ereliability = readProb(graphFile, flowFile, ntests)

function SubproblemH(kth, ntests, dict_dict, c, iter, dual_alpha, dual_beta, dual_pi, dual_eta, sign, quit, init, SPtime)

    # create model
    # model = Model(CPLEX.Optimizer)

    spmodel = Model(Gurobi.Optimizer)
    set_optimizer_attribute(spmodel, "Threads", 1)
    set_optimizer_attribute(spmodel, "OutputFlag", 0)  # 关闭输出

    # def variable
    @variable(spmodel, x[1:nV, 1:nS], Bin)  # x[v, kth, s]
    @variable(spmodel, x0[1:nV], Bin)  # x[v, kth]
    @variable(spmodel, y[1:nV], Bin)  # y[v]
    @variable(spmodel, theta[1:nS+1] >= 0) # theta[kth, s]

    @variable(spmodel, theta_L >= 0) # theta[L, kth]
    @variable(spmodel, theta_N >= 0) # theta[N, kth]
    @variable(spmodel, rel_L)   #
    @variable(spmodel, rel_N)   #

    @variable(spmodel, lambda[1, 1:nS+1, 1:nP] >= 0) # r[kth, s, p]
    @variable(spmodel, z[1:nL, 1:nS+1, 1:nP], Bin)   # z[L, kth, s, p]
    @variable(spmodel, z0[1:nL], Bin)                # z[L, kth]
    @variable(spmodel, r[1:nL, 1:nS+1, 1:nP] >= 0)   # r[L, kth, s, p]
    @variable(spmodel, R[1:nL] >= 0)                 # R[L, kth]
    @variable(spmodel, R0[1:nV] >= 0)                # R[v, kth]

    # add bound constraint
    @constraint(spmodel, r .<= 1)
    @constraint(spmodel, lambda .<= 1)

    VNFplacementforsp(spmodel, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0, kth)

    TrafficroutingMILPforsp(spmodel, x, r, z, V, L, PairC, nV, nI, nL, nK, nS, nP, z0, kth)

    Linkcapforsp(spmodel, r, R, DataRate, Link_cap, nL, nK, kth)

    E2EreliabilityConsforsp(spmodel, x, x0, r, z, z0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, nP, kth)

    E2EdelayConsforsp(spmodel, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP, kth)

    if iter == 1
        @objective(spmodel, Min, sum(y)+0.005sum(R))
    elseif init == 1
        @objective(spmodel, Max, dual_alpha[kth] + sum((dual_beta .- 0.005) .* R) + sum(dual_pi[:, kth] .* x0) + sum(dual_eta .* R0)) 
    else 
        @constraint(spmodel, dual_alpha[kth] + sum((dual_beta .- 0.005) .* R) + sum(dual_pi[:, kth] .* x0) + sum(dual_eta .* R0) >= 1e-3)
        @objective(spmodel, Min, 0)
    end

    optimize!(spmodel)

    status = JuMP.termination_status(spmodel)

    if status == MOI.OPTIMAL
        solve_time = MOI.get(spmodel, MOI.SolveTimeSec())
        SPtime += solve_time
        reducedcost = value.(dual_alpha[kth] + sum((dual_beta .- 0.005) .* R) + sum(dual_pi[:, kth] .* x0) + sum(dual_eta .* R0))
        println("Reduced cost of instance$ntests Subproblem$kth in iter$iter : $reducedcost")
        objective_value = JuMP.objective_value(spmodel)
        #println("Objective Value of instance$ntests Subproblem$kth in iter$iter : $objective_value")
        # if init == 1 
        c[kth] += 1

        solution_dict = Dict{Symbol, Any}()

        solution_dict[:x] = value.(x)
        solution_dict[:x0] = value.(x0)
        solution_dict[:y] = value.(y)
        solution_dict[:theta] = value.(theta)
        solution_dict[:theta_L] = value(theta_L)
        solution_dict[:theta_N] = value(theta_N)
        solution_dict[:rel_L] = value(rel_L)
        solution_dict[:rel_N] = value(rel_N)
        solution_dict[:lambda] = value.(lambda)
        solution_dict[:z] = value.(z)
        solution_dict[:z0] = value.(z0)
        solution_dict[:r] = value.(r)
        solution_dict[:R] = value.(R)
        solution_dict[:R0] = value.(R0)
    
        # println("Values of R: $(solution_dict[:R])")
        # println("Values of x0: $(solution_dict[:x0])")
        # println("Values of R0: $(solution_dict[:R0])")
        #dict_dict = Dict("($kth,$(c[kth]))" => copy(solution_dict))
        # println("key:$kth,$(c[kth])")
        dict_dict[(kth, c[kth])] = copy(solution_dict)
        # 输出字典的键和值
        # println("Keys of the dictionary: $(keys(dict_dict))")
        # println("Values of the dictionary: $(values(dict_dict))")
        # else 
        #     sign[kth] = 1
        #     println("no new patterns in Subproblem$kth in iteraion$iter !")
        # end
    else
        if init == 1
            # 无解
            quit = 1
            println("instance$ntests Subproblem$kth no feasible solution!")
            println("origin problem no feasible solution!")
        else
            sign[kth] = 1
            println("no new patterns in Subproblem$kth in iteraion$iter !")
        end
    end    
    return quit, SPtime 

end

function Subproblem(kth, ntests, dict_dict, c, iter, dual_alpha, dual_beta, dual_pi, dual_eta, sign, quit, init, SPtime)

    # create model
    # model = Model(CPLEX.Optimizer)

    spmodel = Model(Gurobi.Optimizer)
    set_optimizer_attribute(spmodel, "Threads", 1)
    set_optimizer_attribute(spmodel, "OutputFlag", 0)  # 关闭输出

    # def variable
    @variable(spmodel, x[1:nV, 1:nS], Bin)  # x[v, kth, s]
    @variable(spmodel, x0[1:nV], Bin)  # x[v, kth]
    @variable(spmodel, y[1:nV], Bin)  # y[v]
    @variable(spmodel, theta[1:nS+1] >= 0) # theta[kth, s]

    @variable(spmodel, theta_L >= 0) # theta[L, kth]
    @variable(spmodel, theta_N >= 0) # theta[N, kth]
    @variable(spmodel, rel_L)   #
    @variable(spmodel, rel_N)   #

    @variable(spmodel, lambda[1, 1:nS+1, 1:nP] >= 0) # r[kth, s, p]
    @variable(spmodel, z[1:nL, 1:nS+1, 1:nP], Bin)   # z[L, kth, s, p]
    @variable(spmodel, z0[1:nL], Bin)                # z[L, kth]
    @variable(spmodel, r[1:nL, 1:nS+1, 1:nP] >= 0)   # r[L, kth, s, p]
    @variable(spmodel, R[1:nL] >= 0)                 # R[L, kth]
    @variable(spmodel, R0[1:nV] >= 0)                # R[v, kth]

    # add bound constraint
    @constraint(spmodel, r .<= 1)
    @constraint(spmodel, lambda .<= 1)

    VNFplacementforsp(spmodel, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0, kth)

    TrafficroutingMILPforsp(spmodel, x, r, z, V, L, PairC, nV, nI, nL, nK, nS, nP, z0, kth)

    Linkcapforsp(spmodel, r, R, DataRate, Link_cap, nL, nK, kth)

    E2EreliabilityConsforsp(spmodel, x, x0, r, z, z0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, nP, kth)

    E2EdelayConsforsp(spmodel, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP, kth)

    if iter == 1
        @objective(spmodel, Min, sum(y)+0.005sum(R))
    else
        @objective(spmodel, Max, dual_alpha[kth] + sum((dual_beta .- 0.005) .* R) + sum(dual_pi[:, kth] .* x0) + sum(dual_eta .* R0)) 
    end

    optimize!(spmodel)

    status = JuMP.termination_status(spmodel)

    if status == MOI.OPTIMAL
        solve_time = MOI.get(spmodel, MOI.SolveTimeSec())
        SPtime += solve_time

        objective_value = JuMP.objective_value(spmodel)
        println("Objective Value of instance$ntests Subproblem$kth in iter$iter : $objective_value")
        if objective_value > 1e-3 || init == 1 
            c[kth] += 1

            solution_dict = Dict{Symbol, Any}()

            solution_dict[:x] = value.(x)
            solution_dict[:x0] = value.(x0)
            solution_dict[:y] = value.(y)
            solution_dict[:theta] = value.(theta)
            solution_dict[:theta_L] = value(theta_L)
            solution_dict[:theta_N] = value(theta_N)
            solution_dict[:rel_L] = value(rel_L)
            solution_dict[:rel_N] = value(rel_N)
            solution_dict[:lambda] = value.(lambda)
            solution_dict[:z] = value.(z)
            solution_dict[:z0] = value.(z0)
            solution_dict[:r] = value.(r)
            solution_dict[:R] = value.(R)
            solution_dict[:R0] = value.(R0)
        
            # println("Values of R: $(solution_dict[:R])")
            # println("Values of x0: $(solution_dict[:x0])")
            # println("Values of R0: $(solution_dict[:R0])")
            #dict_dict = Dict("($kth,$(c[kth]))" => copy(solution_dict))
            # println("key:$kth,$(c[kth])")
            dict_dict[(kth, c[kth])] = copy(solution_dict)
            # 输出字典的键和值
            # println("Keys of the dictionary: $(keys(dict_dict))")
            # println("Values of the dictionary: $(values(dict_dict))")
        else 
            sign[kth] = 1
            println("no new patterns in Subproblem$kth in iteraion$iter !")
        end
    else
        # 无解
        quit = 1
        println("instance$ntests Subproblem$kth no feasible solution!")
        println("origin problem no feasible solution!")
    end    
    return quit, SPtime 

end

function SubproblemLP(kth, ntests, dict_dict, c, iter, dual_alpha, dual_beta, dual_pi, dual_eta, sign, quit, init, SPtime, condi)

    # create model
    # model = Model(CPLEX.Optimizer)
    
    spmodel = Model(Gurobi.Optimizer)
    set_optimizer_attribute(spmodel, "Threads", 1)
    set_optimizer_attribute(spmodel, "OutputFlag", 0)  # 关闭输出

    # def variable
    @variable(spmodel, x[1:nV, 1:nS])  # x[v, kth, s]
    @variable(spmodel, x0[1:nV])  # x[v, kth]
    @variable(spmodel, y[1:nV])  # y[v]
    @variable(spmodel, theta[1:nS+1] >= 0) # theta[kth, s]

    @variable(spmodel, theta_L >= 0) # theta[L, kth]
    @variable(spmodel, theta_N >= 0) # theta[N, kth]
    @variable(spmodel, rel_L)   #
    @variable(spmodel, rel_N)   #

    @variable(spmodel, lambda[1, 1:nS+1, 1:nP] >= 0) # r[kth, s, p]
    @variable(spmodel, z[1:nL, 1:nS+1, 1:nP])   # z[L, kth, s, p]
    @variable(spmodel, z0[1:nL])                # z[L, kth]
    @variable(spmodel, r[1:nL, 1:nS+1, 1:nP] >= 0)   # r[L, kth, s, p]
    @variable(spmodel, R[1:nL] >= 0)                 # R[L, kth]
    @variable(spmodel, R0[1:nV] >= 0)                # R[v, kth]

    # add bound constraint
    @constraint(spmodel, x .<= 1)
    @constraint(spmodel, x0 .<= 1)
    @constraint(spmodel, y .<= 1)
    @constraint(spmodel, z .<= 1)
    @constraint(spmodel, z0 .<= 1)
    @constraint(spmodel, r .<= 1)
    @constraint(spmodel, lambda .<= 1)

    VNFplacementforsp(spmodel, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0, kth)

    TrafficroutingMILPforsp(spmodel, x, r, z, V, L, PairC, nV, nI, nL, nK, nS, nP, z0, kth)

    Linkcapforsp(spmodel, r, R, DataRate, Link_cap, nL, nK, kth)

    E2EreliabilityConsforsp(spmodel, x, x0, r, z, z0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, nP, kth)

    E2EdelayConsforsp(spmodel, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP, kth)

    if iter == 1
        @objective(spmodel, Min, sum(y)+0.005sum(R))
    else
        @objective(spmodel, Max, dual_alpha[kth] + sum((dual_beta .- 0.005) .* R) + sum(dual_pi[:, kth] .* x0) + sum(dual_eta .* R0)) 
    end

    optimize!(spmodel)

    status = JuMP.termination_status(spmodel)

    if status == MOI.OPTIMAL

        #######把几个整数变量的值存起来######
        sol_x = value.(x)
        sol_x0 = value.(x0)
        sol_y = value.(y)
        sol_z = value.(z)
        sol_z0 = value.(z0)

        solve_time = MOI.get(spmodel, MOI.SolveTimeSec())
        SPtime += solve_time

        objective_value = JuMP.objective_value(spmodel)
        println("Objective Value of instance$ntests Subproblem$kth LP relaxation in iter$iter : $objective_value")
        if objective_value <= 1e-3 && init == 0
            condi = 1
            sign[kth] = 1
            println("no new patterns in Subproblem$kth in iteraion$iter !")
        elseif are_arrays_integers(sol_x, sol_x0, sol_y, sol_z, sol_z0)
            println("Subproblem$kth LP in iteraion$iter return integer solution:")
            println("x=$sol_x, x0=$sol_x0, y=$sol_y, z=$sol_z, z0=$sol_z0")

            condi = 1
            c[kth] += 1

            solution_dict = Dict{Symbol, Any}()

            solution_dict[:x] = value.(x)
            solution_dict[:x0] = value.(x0)
            solution_dict[:y] = value.(y)
            solution_dict[:theta] = value.(theta)
            solution_dict[:theta_L] = value(theta_L)
            solution_dict[:theta_N] = value(theta_N)
            solution_dict[:rel_L] = value(rel_L)
            solution_dict[:rel_N] = value(rel_N)
            solution_dict[:lambda] = value.(lambda)
            solution_dict[:z] = value.(z)
            solution_dict[:z0] = value.(z0)
            solution_dict[:r] = value.(r)
            solution_dict[:R] = value.(R)
            solution_dict[:R0] = value.(R0)
        
            dict_dict[(kth, c[kth])] = copy(solution_dict)

        end
    else
        # 无解
        quit = 1
        println("instance$ntests Subproblem$kth no feasible solution!")
        println("origin problem no feasible solution!")
    end    

    return quit, SPtime, condi 

end

function SubproblemNovelLP(kth, ntests, dict_dict, c, iter, dual_alpha, dual_beta, dual_pi, dual_eta, sign, quit, init, SPLPtime, condi, Novelcolumn, warmstart_dict, model_dict)

    # create model
    # model = Model(CPLEX.Optimizer)
    
    spmodel = Model(Gurobi.Optimizer)
    set_optimizer_attribute(spmodel, "Threads", 1)
    set_optimizer_attribute(spmodel, "OutputFlag", 0)  # 关闭输出

    # def variable
    @variable(spmodel, x[1:nV, 1:nS])  # x[v, kth, s]
    @variable(spmodel, x0[1:nV])  # x[v, kth]
    @variable(spmodel, y[1:nV])  # y[v]
    @variable(spmodel, theta[1:nS+1] >= 0) # theta[kth, s]

    @variable(spmodel, theta_L >= 0) # theta[L, kth]
    @variable(spmodel, theta_N >= 0) # theta[N, kth]
    @variable(spmodel, rel_L)   #
    @variable(spmodel, rel_N)   #

    @variable(spmodel, r[1:nL, 1:nS+1] >= 0)   # r[L, kth, s, 1]
    @variable(spmodel, r0[1:nL] >= 0)
    @variable(spmodel, R[1:nL] >= 0)                 # R[L, kth]
    @variable(spmodel, R0[1:nV] >= 0)                # R[v, kth]

    # add bound constraint
    @constraint(spmodel, x .<= 1)
    @constraint(spmodel, x0 .<= 1)
    @constraint(spmodel, y .<= 1)
    #@constraint(spmodel, z .<= 1)
    #@constraint(spmodel, z0 .<= 1)
    @constraint(spmodel, r .<= 1)
    @constraint(spmodel, r0 .<= 1)
    #@constraint(spmodel, lambda .<= 1)

    VNFplacementforspNovelLP(spmodel, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0, kth)

    TrafficroutingMILPforspNovelLP(spmodel, x, r, V, L, PairC, nV, nI, nL, nK, nS, kth)

    LinkcapforspNovelLP(spmodel, r, R, DataRate, Link_cap, nL, nK, kth)

    E2EreliabilityConsforspNovelLP(spmodel, x, x0, r, r0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, kth)

    E2EdelayConsforspNovelLP(spmodel, x, r, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, kth)

    if iter == 1
        @objective(spmodel, Min, sum(y)+0.005sum(R))
    elseif iter > 2
        # set_start_value.(x, warmstart_dict[kth][:x])
        # set_start_value.(x0, warmstart_dict[kth][:x0])
        # set_start_value.(y, warmstart_dict[kth][:y])
        # set_start_value.(r, warmstart_dict[kth][:r])
        # set_start_value.(r0, warmstart_dict[kth][:r0])
        # set_start_value.(theta, warmstart_dict[kth][:theta])
        # set_start_value.(theta_L, warmstart_dict[kth][:theta_L])
        # set_start_value.(theta_N, warmstart_dict[kth][:theta_N])
        # set_start_value.(rel_L, warmstart_dict[kth][:rel_L])
        # set_start_value.(rel_N, warmstart_dict[kth][:rel_N])
        # set_start_value.(R, warmstart_dict[kth][:R])
        # set_start_value.(R0, warmstart_dict[kth][:R0])

        @objective(spmodel, Max, dual_alpha[kth] + sum((dual_beta .- 0.005) .* R) + sum(dual_pi[:, kth] .* x0) + sum(dual_eta .* R0)) 
    else
        @objective(spmodel, Max, dual_alpha[kth] + sum((dual_beta .- 0.005) .* R) + sum(dual_pi[:, kth] .* x0) + sum(dual_eta .* R0)) 
    end

    optimize!(spmodel)

    status = JuMP.termination_status(spmodel)

    if status == MOI.OPTIMAL

        sol_x = value.(x)
        sol_x0 = value.(x0)
        sol_y = value.(y)
        sol_r = value.(r)
        sol_r0 = value.(r0)

        start_dict = Dict{Symbol, Any}()
        start_dict[:x] = value.(x)
        start_dict[:x0] = value.(x0)
        start_dict[:y] = value.(y)
        start_dict[:theta] = value.(theta)
        start_dict[:theta_L] = value(theta_L)
        start_dict[:theta_N] = value(theta_N)
        start_dict[:rel_L] = value(rel_L)
        start_dict[:rel_N] = value(rel_N)
        start_dict[:r] = value.(r)
        start_dict[:r0] = value.(r0)
        start_dict[:R] = value.(R)
        start_dict[:R0] = value.(R0)

        warmstart_dict[kth] = copy(start_dict)

        solve_time = MOI.get(spmodel, MOI.SolveTimeSec())
        SPLPtime += solve_time

        objective_value = JuMP.objective_value(spmodel)
        println("Objective Value of instance$ntests Subproblem$kth LP relaxation in iter$iter : $objective_value")
        if objective_value <= 1e-3 && init == 0
            condi = 1
            sign[kth] = 1
            println("no new patterns in Subproblem$kth in iteraion$iter !")
        elseif are_arrays_integers(sol_x, sol_x0, sol_y, sol_r, sol_r0)
            println("Subproblem$kth NolvelLP in iteraion$iter return integer solution:")
           # println("x=$sol_x, x0=$sol_x0, y=$sol_y, r=$sol_r, r0=$sol_r0")
            Novelcolumn += 1
            condi = 1
            c[kth] += 1

            solution_dict = Dict{Symbol, Any}()

            solution_dict[:x] = value.(x)
            solution_dict[:x0] = value.(x0)
            solution_dict[:y] = value.(y)
            solution_dict[:theta] = value.(theta)
            solution_dict[:theta_L] = value(theta_L)
            solution_dict[:theta_N] = value(theta_N)
            solution_dict[:rel_L] = value(rel_L)
            solution_dict[:rel_N] = value(rel_N)
            #solution_dict[:lambda] = value.(lambda)
            #solution_dict[:z] = value.(z)
            #solution_dict[:z0] = value.(z0)
            solution_dict[:r] = value.(r)
            solution_dict[:r0] = value.(r0)
            solution_dict[:R] = value.(R)
            solution_dict[:R0] = value.(R0)
        
            dict_dict[(kth, c[kth])] = copy(solution_dict)

        end
    else
        # 无解
        quit = 1
        println("instance$ntests Subproblem$kth no feasible solution!")
        println("origin problem no feasible solution!")
    end    

    return quit, SPtime, condi, Novelcolumn

end

function SubproblemNovelLPWarmstart(kth, ntests, dict_dict, c, iter, dual_alpha, dual_beta, dual_pi, dual_eta, sign, quit, init, SPLPtime, condi, Novelcolumn, warmstart_dict, model_dict)

    # create model
    # model = Model(CPLEX.Optimizer)
    if iter == 2
        spmodel = Model(Gurobi.Optimizer)
        
        set_optimizer_attribute(spmodel, "Threads", 1)
        set_optimizer_attribute(spmodel, "OutputFlag", 0)  # 关闭输出
        set_optimizer_attribute(spmodel, "MIPGap", 0)

        # def variable
        @variable(spmodel, x[1:nV, 1:nS])  # x[v, kth, s]
        @variable(spmodel, x0[1:nV])  # x[v, kth]
        @variable(spmodel, y[1:nV])  # y[v]
        @variable(spmodel, theta[1:nS+1] >= 0) # theta[kth, s]

        @variable(spmodel, theta_L >= 0) # theta[L, kth]
        @variable(spmodel, theta_N >= 0) # theta[N, kth]
        @variable(spmodel, rel_L)   #
        @variable(spmodel, rel_N)   #

        @variable(spmodel, r[1:nL, 1:nS+1] >= 0)   # r[L, kth, s, 1]
        @variable(spmodel, r0[1:nL] >= 0)
        @variable(spmodel, R[1:nL] >= 0)                 # R[L, kth]
        @variable(spmodel, R0[1:nV] >= 0)                # R[v, kth]

        # add bound constraint
        @constraint(spmodel, x .<= 1)
        @constraint(spmodel, x0 .<= 1)
        @constraint(spmodel, y .<= 1)
        #@constraint(spmodel, z .<= 1)
        #@constraint(spmodel, z0 .<= 1)
        @constraint(spmodel, r .<= 1)
        @constraint(spmodel, r0 .<= 1)
        #@constraint(spmodel, lambda .<= 1)

        VNFplacementforspNovelLP(spmodel, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0, kth)

        TrafficroutingMILPforspNovelLP(spmodel, x, r, V, L, PairC, nV, nI, nL, nK, nS, kth)

        LinkcapforspNovelLP(spmodel, r, R, DataRate, Link_cap, nL, nK, kth)

        E2EreliabilityConsforspNovelLP(spmodel, x, x0, r, r0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, kth)

        E2EdelayConsforspNovelLP(spmodel, x, r, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, kth)

        @objective(spmodel, Max, dual_alpha[kth] + sum((dual_beta .- 0.005) .* R) + sum(dual_pi[:, kth] .* x0) + sum(dual_eta .* R0)) 
        
        model_dict[kth] = spmodel

        optimize!(spmodel)

        status = JuMP.termination_status(spmodel)
    elseif iter > 2
        # set_start_value.(model_dict[kth][:x], warmstart_dict[kth][:x])
        # set_start_value.(model_dict[kth][:x0], warmstart_dict[kth][:x0])
        # set_start_value.(model_dict[kth][:y], warmstart_dict[kth][:y])
        # set_start_value.(model_dict[kth][:r], warmstart_dict[kth][:r])
        # set_start_value.(model_dict[kth][:r0], warmstart_dict[kth][:r0])
        # set_start_value.(model_dict[kth][:theta], warmstart_dict[kth][:theta])
        # set_start_value.(model_dict[kth][:theta_L], warmstart_dict[kth][:theta_L])
        # set_start_value.(model_dict[kth][:theta_N], warmstart_dict[kth][:theta_N])
        # set_start_value.(model_dict[kth][:rel_L], warmstart_dict[kth][:rel_L])
        # set_start_value.(model_dict[kth][:rel_N], warmstart_dict[kth][:rel_N])
        # set_start_value.(model_dict[kth][:R], warmstart_dict[kth][:R])
        # set_start_value.(model_dict[kth][:R0], warmstart_dict[kth][:R0])

        @objective(model_dict[kth], Max, dual_alpha[kth] + sum((dual_beta .- 0.005) .* model_dict[kth][:R]) + sum(dual_pi[:, kth] .* model_dict[kth][:x0]) + sum(dual_eta .* model_dict[kth][:R0])) 
        
        optimize!(model_dict[kth])
        
        status = JuMP.termination_status(model_dict[kth])

    end

    if status == MOI.OPTIMAL

        sol_x = value.(model_dict[kth][:x])
        sol_x0 = value.(model_dict[kth][:x0])
        sol_y = value.(model_dict[kth][:y])
        sol_r = value.(model_dict[kth][:r])
        sol_r0 = value.(model_dict[kth][:r0])

        start_dict = Dict{Symbol, Any}()
        start_dict[:x] = value.(model_dict[kth][:x])
        start_dict[:x0] = value.(model_dict[kth][:x0])
        start_dict[:y] = value.(model_dict[kth][:y])
        start_dict[:theta] = value.(model_dict[kth][:theta])
        start_dict[:theta_L] = value(model_dict[kth][:theta_L])
        start_dict[:theta_N] = value(model_dict[kth][:theta_N])
        start_dict[:rel_L] = value(model_dict[kth][:rel_L])
        start_dict[:rel_N] = value(model_dict[kth][:rel_N])
        start_dict[:r] = value.(model_dict[kth][:r])
        start_dict[:r0] = value.(model_dict[kth][:r0])
        start_dict[:R] = value.(model_dict[kth][:R])
        start_dict[:R0] = value.(model_dict[kth][:R0])

        warmstart_dict[kth] = copy(start_dict)

        if iter == 2
            solve_time = MOI.get(spmodel, MOI.SolveTimeSec())
            objective_value = JuMP.objective_value(spmodel)
        else
            solve_time = MOI.get(model_dict[kth], MOI.SolveTimeSec())
            objective_value = JuMP.objective_value(model_dict[kth])
        end
        SPLPtime += solve_time

        
        println("Objective Value of instance$ntests Subproblem$kth LP relaxation in iter$iter : $objective_value")
        if objective_value <= 1e-3 && init == 0
            condi = 1
            sign[kth] = 1
            println("no new patterns in Subproblem$kth in iteraion$iter !")
        elseif are_arrays_integers(sol_x, sol_x0, sol_y, sol_r, sol_r0)
            println("Subproblem$kth NolvelLP in iteraion$iter return integer solution:")
           # println("x=$sol_x, x0=$sol_x0, y=$sol_y, r=$sol_r, r0=$sol_r0")
            Novelcolumn += 1
            condi = 1
            c[kth] += 1

            solution_dict = Dict{Symbol, Any}()

            if iter == 2
                solution_dict[:x] = value.(x)
                solution_dict[:x0] = value.(x0)
                solution_dict[:y] = value.(y)
                solution_dict[:theta] = value.(theta)
                solution_dict[:theta_L] = value(theta_L)
                solution_dict[:theta_N] = value(theta_N)
                solution_dict[:rel_L] = value(rel_L)
                solution_dict[:rel_N] = value(rel_N)
                solution_dict[:r] = value.(r)
                solution_dict[:r0] = value.(r0)
                solution_dict[:R] = value.(R)
                solution_dict[:R0] = value.(R0)
            else
                solution_dict[:x] = value.(model_dict[kth][:x])
                solution_dict[:x0] = value.(model_dict[kth][:x0])
                solution_dict[:y] = value.(model_dict[kth][:y])
                solution_dict[:theta] = value.(model_dict[kth][:theta])
                solution_dict[:theta_L] = value(model_dict[kth][:theta_L])
                solution_dict[:theta_N] = value(model_dict[kth][:theta_N])
                solution_dict[:rel_L] = value(model_dict[kth][:rel_L])
                solution_dict[:rel_N] = value(model_dict[kth][:rel_N])
                solution_dict[:r] = value.(model_dict[kth][:r])
                solution_dict[:r0] = value.(model_dict[kth][:r0])
                solution_dict[:R] = value.(model_dict[kth][:R])
                solution_dict[:R0] = value.(model_dict[kth][:R0])
            end
        
            dict_dict[(kth, c[kth])] = copy(solution_dict)

        end
    else
        # 无解
        quit = 1
        println("instance$ntests Subproblem$kth no feasible solution!")
        println("origin problem no feasible solution!")
    end    

    return quit, SPLPtime, condi, Novelcolumn

end

