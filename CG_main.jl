using JuMP
using CPLEX
using Gurobi
using Printf
using DelimitedFiles
using Statistics
import MathOptInterface
const MOI = MathOptInterface
using Random
using Dates

include("ReadData.jl")
include("constraint.jl")
include("subproblem.jl")
include("masterproblem.jl")

Profun = 3
nP = 2

function CG_main(ntest, K, capacityparam)

    file_path = "/share/home/sunqi/GitPublish/ns-algorithm/julia_code/log/CGWS-$K-$capacityparam-$ntest.txt"
    global file = open(file_path, "a")
    
    global SPtime = 0
    global SPLPtime = 0
    global MPtime = 0    

    global ntests = ntest
    graphFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-$K-$capacityparam-graph_info-$ntests"
    flowFile = "/share/home/sunqi/GitPublish/ns-data/instances/fish/fish-$K-$capacityparam-flow_info-$ntests"

    global nI, nL, nV, nF, L, V, Link_cap, Link_reliability, Link_delay, Node_cap, Node_delay, Node_reliability, outgoing_cap, incoming_cap, V_processfunction, nK, nS, PairC, DataRate, FunctionChainC, E2EDelay, E2Ereliability = readProb(graphFile, flowFile, ntests)


    global dict_dict = Dict{Tuple{Int, Int}, Any}()
    global warmstart_dict = Dict{Int, Any}()
    global model_dict = Dict{Int, Any}()
    global dual_alpha = zeros(nK)
    global dual_beta = zeros(nL)
    global dual_pi = zeros(nV, nK)
    global dual_eta = zeros(nV)
    global Novelcolumn = 0
    global c = zeros(Int, nK)
    global init = 1         # 用于判断松弛主问题是否已经找到解
    global quit = 0         # 用于判断是否存在一个子问题无解

    # 获取开始时间
    start_time = now()

    for iter = 1 : 20
        
        sign = zeros(nK)

        if iter != 1
            # 先解子问题的松弛
            for kth = 1 : nK
                condi = 0
                quit, SPLPtime, condi, Novelcolumn = SubproblemNovelLPWarmstart(kth, ntests, dict_dict, c, iter, dual_alpha, dual_beta, dual_pi, dual_eta, sign, quit, init, SPLPtime, condi, Novelcolumn, warmstart_dict, model_dict)
                if condi == 0
                    quit, SPtime = Subproblem(kth, ntests, dict_dict, c, iter, dual_alpha, dual_beta, dual_pi, dual_eta, sign, quit, init, SPtime)
                end
            end
        else
            # 子问题-生成新列
            for kth = 1 : nK
                quit, SPtime = Subproblem(kth, ntests, dict_dict, c, iter, dual_alpha, dual_beta, dual_pi, dual_eta, sign, quit, init, SPtime)
            end
        end


        if quit == 1
            println(file,"CG-$K-$capacityparam instance$ntests has no feasible solution!")
            break
        end

        if sum(sign) == nK || iter == 20
            if sum(sign) == nK
                println("CG instance $ntests Complete!")
            else
                println(file, "CG-$K-$capacityparam fail to complete instance $ntests in $iter iterations")
            end
            masterproblemMIP(ntests)
            println(file, "CG instance$ntests finished at the $iter-th iteraion ")
            break
        end

        # 求解松弛主问题
        global dual_alpha, dual_beta, dual_pi, dual_eta, init, MPtime = masterproblem(iter, ntests, init, MPtime)

    end

    # 获取结束时间
    end_time = now()

    # 计算时间差
    elapsed_time = end_time - start_time

    #elapsed_seconds = value(elapsed_time)

    println(file, "CG instances $ntests finished in $elapsed_time")
    println(file, "CG instance $ntests SP time is $SPtime")
    println(file, "CG instance $ntests MILP time is $SPLPtime")
    println(file, "CG instance $ntests MP time is $MPtime")
    println(file, "CG instance $ntests column number is $(sum(c))")
    println(file, "CG instance $ntests Novelinteger number is $Novelcolumn")
    
end

# 从命令行获取参数
args = parse.(Int, ARGS)

# 调用 main 函数

CG_main(args...)


