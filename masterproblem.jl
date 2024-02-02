function masterproblem(iter, ntests, init, MPtime)
 
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "Threads", 1)
    if init == 1
        set_optimizer_attribute(model, "Presolve", 0)
    else
        set_optimizer_attribute(model, "Presolve", 1)
    end
    set_optimizer_attribute(model, "OutputFlag", 0)  # 关闭输出

    nkth = 1 : nK

    # create variable
    @variable(model, t[i=1:length(nkth), j=1:c[nkth[i]]]>=0)
    @variable(model, y[1:nV]>=0)

    # add constraints
    @constraint(model, y .<= 1)

    consalpha, consbeta, conspi, conseta, Lhs1, model = ConsforMP(model, t, y, c, dict_dict, Link_cap, Node_cap, nK, nL)

    # println("?????????????????????????????/")
    # println(Lhs1)

    # Objective function
    @objective(model, Min, sum(y)+0.005*sum(Lhs1))

    # solve RLM
    optimize!(model)

    status = JuMP.termination_status(model)
    dualstatus = dual_status(model)

    if status == MOI.OPTIMAL

        init = 0
        #println(file, "CG instance$ntests Masterproblem LP finds feasible solution at the $iter th iteraion")
        solve_time = MOI.get(model, MOI.SolveTimeSec())
        MPtime += solve_time
        objective_valueM = JuMP.objective_value(model)
        println("Objective Value of instance$ntests Masterproblem in iter$iter : $objective_valueM")

        # get dual
        dual_alpha = [dual(consalpha[k]) for k in 1:nK]
        dual_beta = [dual(consbeta[l]) for l in 1:nL]
        dual_pi = [dual(conspi[v, k]) for v in 1:nV, k in 1:nK]
        dual_eta = [dual(conseta[v]) for v in 1:nV]

        # print result
        # println("Optimal solution of t in iter$iter:")
        # println(value.(t))
        # println("Optimal solution of y in iter$iter:")
        # println(value.(y))
        # println("Dual values for constraints alpha:")
        # println(dual_alpha)
        # println("Dual values for constraints beta:")
        # println(dual_beta)
        # println("Dual values for constraints pi:")
        # println(dual_pi)
        # println("Dual values for constraints eta:")
        # println(dual_eta)
        # println(size(dual_eta))

    else

        solve_time = MOI.get(model, MOI.SolveTimeSec())
        MPtime += solve_time
        println("Master problem has no feasible solution in iteraion$iter yet! ")
        println("status is $status")
        println("dual status is $dualstatus")
        # get dual rays
        dual_alpha = [dual(consalpha[k]) for k in 1:nK]
        dual_beta = [dual(consbeta[l]) for l in 1:nL]
        dual_pi = [dual(conspi[v, k]) for v in 1:nV, k in 1:nK]
        dual_eta = [dual(conseta[v]) for v in 1:nV]

        # println("Dual values for constraints alpha:")
        # println(dual_alpha)
        # println("Dual values for constraints beta:")
        # println(dual_beta)
        # println("Dual values for constraints pi:")
        # println(dual_pi)
        # println("Dual values for constraints eta:")
        # println(dual_eta)
    end
    return dual_alpha, dual_beta, dual_pi, dual_eta, init, MPtime

end

function masterproblemMIP(ntests)
 
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "Threads", 1)
    set_optimizer_attribute(model, "Presolve", 1)
    set_optimizer_attribute(model, "OutputFlag", 0)  # 关闭输出

    nkth = 1 : nK

    # create variable
    @variable(model, t[i=1:length(nkth), j=1:c[nkth[i]]]>=0, Bin)
    @variable(model, y[1:nV]>=0, Bin)

    # add constraints

    consalpha, consbeta, conspi, conseta, Lhs1, model = ConsforMP(model, t, y, c, dict_dict, Link_cap, Node_cap, nK, nL)

    # println("?????????????????????????????/")
    # println(Lhs1)

    # Objective function
    @objective(model, Min, sum(y)+0.005*sum(Lhs1))

    # solve RLM
    optimize!(model)

    status = JuMP.termination_status(model)

    if status == MOI.OPTIMAL

        solve_time = MOI.get(model, MOI.SolveTimeSec())

        objective_valueF = JuMP.objective_value(model)
        println(file, "CG objective Value of instance$ntests : $objective_valueF")
        println(file, "CG instance$ntests MPMIP solve time is $solve_time")

        # print result
        # println("Optimal solution of t:")
        # println(value.(t))
        # println("Optimal solution of y:")
        # println(value.(y))
        # println("Dual values for constraints alpha:")
        # println(dual_alpha)
        # println("Dual values for constraints beta:")
        # println(dual_beta)
        # println("Dual values for constraints pi:")
        # println(dual_pi)
        # println("Dual values for constraints eta:")
        # println(dual_eta)
        #println(size(dual_eta))
    else
        println(file, "CG failed to find feasible solution of instance $ntests")
    end

end


