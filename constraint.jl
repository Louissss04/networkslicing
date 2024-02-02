function VNFplacement(model, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0)

    ## this function add constraints (1)-(4) in paper

    # no provide constraint (which doesn't appear in the paper) 
    for k = 1:nK
        for s = 1:nS
            for v = 2:nV      # here we always assume that node V[1] can process all service functions
                curprocessfunctions = V_processfunction[1:Profun, v]
                if sum(FunctionChainC[s, k] .== curprocessfunctions) == 0
                    @constraint(model, x[v, k, s] == 0)
                end
            end
        end
    end

    # every function should be processed (and only once)
    for k = 1:nK
        for s = 1:nS
            addition = 0
            for v = 1:nV
                addition += x[v, k, s]
            end
            @constraint(model, addition == 1)
        end
    end

    # x[v, k, s] and x[v, k]
    for k = 1:nK
        for s = 1:nS
            for v = 1:nV
                @constraint(model, x[v, k, s] <= x0[v, k])
            end
        end
    end

    # x[v,k] and y[v]
    for k = 1:nK
        for v = 1:nV
            @constraint(model, x0[v, k] <= y[v])
        end
    end

    # node capacity constraint
    for v = 1:nV
        for k = 1:nK
            addition = 0
            for s = 1:nS
                addition += DataRate[k] * x[v, k, s]
            end
            @constraint(model, addition == R0[v, k])
        end
        @constraint(model, sum(R0[v, :]) <= Node_cap[v] * y[v])
    end

    # aggregate 1 and 4
    addition = 0
    addition2 = 0
    for v = 1:nV
        addition += Node_cap[v] * y[v]
    end

    for k = 1:nK
        for s = 1:nS
            addition2 += DataRate[k]
        end
    end

    @constraint(model, addition2 <= addition)

        # node capacity constraint
        # for v = 1:nV
        #     addition = 0
        #     for k = 1:nK
        #         for s = 1:nS
        #             addition += DataRate[k] * x[v, k, s]
        #         end
        #     end
        #     @constraint(model, addition <= Node_cap[v] * y[v])
        # end

    return model 
end

function TrafficroutingMINLP(model, x, r, z, lambda, V, L, PairC, nV, nI, nL, nK, nS, nP, z0)
    
    ## this function add constraints minlp in paper  

    ## Data rate constraints
    for k = 1:nK
        for s = 1:(nS+1)
            addition = 0
            for p = 1:nP
                addition += lambda[1, k, s, p]
            end
            @constraint(model, addition == 1)
        end
    end

    ## r z relation constraint
    for k = 1:nK
        for l = 1:nL
            for s = 1:(nS+1)
                for p = 1:nP
                    @constraint(model, r[l, k, s, p] <= z[l, k, s, p])
                    @constraint(model, r[l, k, s, p] >= z[l, k, s, p] + lambda[1, k, s, p] - 1)
                end
            end
        end
    end

    ## SFC constraints
    for k = 1:nK
        for s = 1:(nS+1)
            for i = 1:nI
                outgoinglinks = findall(L[:, 1] .== i)
                incominglinks = findall(L[:, 2] .== i)
                # addition2 = zeros(1, 1, 1, nP)
                # addition4 = zeros(1, 1, 1, nP)

                # for t = 1:length(outgoinglinks)
                #     addition2 .= addition2 .+ z[outgoinglinks[t], k, s, :]
                # end
                
                # for t = 1:length(incominglinks)
                #     addition4 .= addition4 .+ z[incominglinks[t], k, s, :]
                # end

                addition2 = sum(z[outgoinglinks, k, s, :], dims=1)
                addition4 = sum(z[incominglinks, k, s, :], dims=1)

                # println(size(addition2))
                
                v = sum((1:nV) .* (i .== V))
                
                # @constraint(model, addition2 .<= ones(1, 1, 1, nP))
                # @constraint(model, addition2 .<= 1)
                
                if s == 1
                    if i == PairC[k, 1]
                        # @constraint(model, addition2 .== ones(1, 1, 1, nP))
                        # @constraint(model, addition4 .== zeros(1, 1, 1, nP))
                        @constraint(model, addition2 .== 1)     # 情况1
                        @constraint(model, addition4 .== 0)
                    elseif v > 0
                        # @constraint(model, addition4 .- addition2 .== ones(1, 1, 1, nP) .* x[v, k, s])
                        @constraint(model, addition4 .- addition2 .== x[v, k, s])    # 情况2
                    else
                        @constraint(model, addition2 .== addition4)    # 情况6
                    end
                elseif s == (nS+1)
                    if i == PairC[k, 2]
                        # @constraint(model, addition2 .== zeros(1, 1, 1, nP))
                        # @constraint(model, addition4 .== ones(1, 1, 1, nP))
                        @constraint(model, addition2 .== 0)                       # 情况5
                        @constraint(model, addition4 .== 1)
                    elseif v > 0
                        #@constraint(model, addition4 .- addition2 .== ones(1, 1, 1, nP) .* (-x[v, k, s-1]))
                        @constraint(model, addition4 .- addition2 .== (-x[v, k, s-1]))          # 情况4
                    else
                        @constraint(model, addition2 .== addition4)                # 情况6
                    end
                else
                    if v > 0
                        #@constraint(model, addition4 .- addition2 .== ones(1, 1, 1, nP) .* (x[v, k, s] - x[v, k, s-1]))
                        @constraint(model, addition4 .- addition2 .== (x[v, k, s] - x[v, k, s-1]))       # 情况3
                    else
                        @constraint(model, addition2 .== addition4)            # 情况6
                    end
                end
            end
        end
    end

        # Valid inequalities
        # for s = 1:(nS + 1)
        #     for l = 1:nL
        #         for k = 1:nK
        #             addition = 0
        #             for p =1:nP
        #                 addition = addition + r[l, k, s, p]
        #             end
        #             @constraint(model, addition <= z0[l, k])
        #         end
        #     end
            # for p = 1:nP
            #     #@constraint(model, z[:, :, s, p] .<= z0[:, :])
            #     addition .= addition .+ r[:, :, s, p]
            # end
            #@constraint(model, addition .<= z0[:, :])
        # end

    return model
end

function TrafficroutingMILP(model, x, r, z, V, L, PairC, nV, nI, nL, nK, nS, nP, z0)

    ## SFC constraints
    for k = 1:nK
        for s = 1:(nS+1)
            for i = 1:nI
                outgoinglinks = findall(L[:, 1] .== i)
                incominglinks = findall(L[:, 2] .== i)
                # addition2 = zeros(1, 1, 1, nP)
                # addition4 = zeros(1, 1, 1, nP)

                # for t = 1:length(outgoinglinks)
                #     addition2 .= addition2 .+ z[outgoinglinks[t], k, s, :]
                # end
                
                # for t = 1:length(incominglinks)
                #     addition4 .= addition4 .+ z[incominglinks[t], k, s, :]
                # end

                addition2 = sum(z[outgoinglinks, k, s, :], dims=1)
                addition4 = sum(z[incominglinks, k, s, :], dims=1)

                # println(size(addition2))
                
                v = sum((1:nV) .* (i .== V))
                
                # @constraint(model, addition2 .<= ones(1, 1, 1, nP))
                @constraint(model, addition2 .<= 1)
                
                if s == 1
                    if i == PairC[k, 1]
                        # @constraint(model, addition2 .== ones(1, 1, 1, nP))
                        # @constraint(model, addition4 .== zeros(1, 1, 1, nP))
                        @constraint(model, addition2 .== 1)     # 情况1
                        @constraint(model, addition4 .== 0)
                    elseif v > 0
                        # @constraint(model, addition4 .- addition2 .== ones(1, 1, 1, nP) .* x[v, k, s])
                        @constraint(model, addition4 .- addition2 .== x[v, k, s])    # 情况2
                    else
                        @constraint(model, addition2 .== addition4)    # 情况6
                    end
                elseif s == (nS+1)
                    if i == PairC[k, 2]
                        # @constraint(model, addition2 .== zeros(1, 1, 1, nP))
                        # @constraint(model, addition4 .== ones(1, 1, 1, nP))
                        @constraint(model, addition2 .== 0)                       # 情况5
                        @constraint(model, addition4 .== 1)
                    elseif v > 0
                        #@constraint(model, addition4 .- addition2 .== ones(1, 1, 1, nP) .* (-x[v, k, s-1]))
                        @constraint(model, addition4 .- addition2 .== (-x[v, k, s-1]))          # 情况4
                    else
                        @constraint(model, addition2 .== addition4)                # 情况6
                    end
                else
                    if v > 0
                        #@constraint(model, addition4 .- addition2 .== ones(1, 1, 1, nP) .* (x[v, k, s] - x[v, k, s-1]))
                        @constraint(model, addition4 .- addition2 .== (x[v, k, s] - x[v, k, s-1]))       # 情况3
                    else
                        @constraint(model, addition2 .== addition4)            # 情况6
                    end
                end
            end
        end
    end

    # rz relation constraint
    for k = 1:nK
        @constraint(model, r[:, k, :, :] .<= z[:, k, :, :])
    end

    for k = 1:nK
        for s = 1:(nS + 1)
            for i = 1:nI
                outgoinglinks = findall(L[:, 1] .== i)
                incominglinks = findall(L[:, 2] .== i)
                addition1 = 0
                addition2 = zeros(1, 1, 1, nP)
                addition3 = 0
                addition4 = zeros(1, 1, 1, nP)
                v = sum((1:nV) .* (i .== V))

                for t = 1:length(outgoinglinks)
                    for p = 1:nP
                        addition1 += r[outgoinglinks[t], k, s, p]
                    end
                    addition2 = addition2 .+ r[outgoinglinks[t], k, s, :]
                end

                for t = 1:length(incominglinks)
                    for p = 1:nP
                        addition3 += r[incominglinks[t], k, s, p]
                    end
                    addition4 = addition4 .+ r[incominglinks[t], k, s, :]

                end

                if s == 1
                    if i == PairC[k, 1]
                        @constraint(model, addition3 - addition1 == -1)      # 15 情况 1
                    elseif v > 0
                        @constraint(model, addition3 - addition1 == x[v, k, s])   # 15 情况 2
                        @constraint(model, addition4 .- addition2 .<= x[v, k, s])   # 17
                    else
                        @constraint(model, addition2 .== addition4)               # 16
                    end
                elseif s == (nS + 1)
                    if i == PairC[k, 2]
                        @constraint(model, addition3 - addition1 == 1)    # 15 情况 5
                    elseif v > 0
                        @constraint(model, addition3 - addition1 == -x[v, k, s - 1])     # 15 情况 4
                        @constraint(model, addition4 .- addition2 .>= -x[v, k, s - 1])    # 18   
                    else
                        @constraint(model, addition2 .== addition4)         # 16
                    end
                else
                    if v > 0
                        @constraint(model, addition3 - addition1 == (x[v, k, s] - x[v, k, s - 1]))         # 15 情况 3   
                        @constraint(model, addition4 .- addition2 .<= x[v, k, s])       # 17
                        @constraint(model, addition4 .- addition2 .>= -x[v, k, s - 1])      # 18
                    else
                        @constraint(model, addition2 .== addition4)         # 16
                    end
                end
            end
        end
    end

    # Valid inequalities
        for s = 1:(nS + 1)
            for l = 1:nL
                for k = 1:nK
                    addition = 0
                    for p =1:nP
                        addition = addition + r[l, k, s, p]
                    end
                    @constraint(model, addition <= z0[l, k])
                end
            end
            # for p = 1:nP
            #     #@constraint(model, z[:, :, s, p] .<= z0[:, :])
            #     addition .= addition .+ r[:, :, s, p]
            # end
            #@constraint(model, addition .<= z0[:, :])
        end

    return model
end


function Linkcap(model, r, R, DataRate, Link_cap, nL, nK)
    # r R related constraint
    for l = 1:nL
        for k = 1:nK
            @constraint(model, R[l, k] == sum(sum(DataRate[k] * r[l, k, :, :])))
        end
    end

    # Link Capacity Constraint
    for l = 1:nL
        @constraint(model, sum(R[l, :]) <= Link_cap[l])
    end

    return model
end

function E2EreliabilityCons(model, x, x0, r, z, z0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, nP)
    # z z0 relation
    for s = 1:(nS + 1)
        #addition = zeros(nL, nK)
        for p = 1:nP
            @constraint(model, z[:, :, s, p] .<= z0[:, :])
           # addition .= addition .+ r[:, :, s, p]
        end
        # @constraint(model, addition .<= z0[:, :])
    end

    for k = 1:nK
        addition = 0
        for l = 1:nL
            addition = addition + log(Link_reliability[l]) * z0[l, k]
        end
        @constraint(model, addition .== rel_L[k])
        addition = 0
        for v = 1:nV
            addition = addition + log(Node_reliability[v]) * x0[v, k]
        end
        @constraint(model, addition .== rel_N[k])

        @constraint(model, rel_L[k] + rel_N[k] >= log(E2Ereliability[k]))
    end

    return model
end

function E2EdelayCons(model, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP)

    for k = 1:nK
        for s = 1:(nS + 1)
            addition2 = 0
            for p = 1:nP
                addition = 0
                for l = 1:nL
                    addition += Link_delay[l] * z[l, k, s, p]
                    addition2 += Link_delay[l] * r[l, k, s, p]
                end
                @constraint(model, addition <= theta[k, s])
            end
            @constraint(model, addition2 <= theta[k, s])
        end
    end

    for k = 1:nK
        @constraint(model, sum(theta[k, :]) == theta_L[k])
    end

    for k = 1:nK
        addition = 0
        for v = 1:nV
            addition += sum(x[v, k, :] .* Node_delay[v])
        end
        @constraint(model, addition == theta_N[k])
    end

    for k = 1:nK
        @constraint(model, theta_L[k] + theta_N[k] <= E2EDelay[k])
    end

    return model
end

function E2EdelayConsforMINLP(model, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP)

    for k = 1:nK
        for s = 1:(nS + 1)
            #addition2 = 0
            for p = 1:nP
                addition = 0
                for l = 1:nL
                    addition += Link_delay[l] * z[l, k, s, p]
                    #addition2 += Link_delay[l] * r[l, k, s, p]
                end
                @constraint(model, addition <= theta[k, s])
            end
           # @constraint(model, addition2 <= theta[k, s])
        end
    end

    for k = 1:nK
        @constraint(model, sum(theta[k, :]) == theta_L[k])
    end

    for k = 1:nK
        addition = 0
        for v = 1:nV
            addition += sum(x[v, k, :] .* Node_delay[v])
        end
        @constraint(model, addition == theta_N[k])
    end

    for k = 1:nK
        @constraint(model, theta_L[k] + theta_N[k] <= E2EDelay[k])
    end

    return model
end

#############以下函数用于列生成算法子问题，即k=1.

function VNFplacementforsp(model, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0, kth)

    ## this function add constraints (1)-(4) in paper

    # no provide constraint (which doesn't appear in the paper) 
    for s = 1:nS
        for v = 2:nV      # here we always assume that node V[1] can process all service functions
            curprocessfunctions = V_processfunction[1:Profun, v]
            if sum(FunctionChainC[s, kth] .== curprocessfunctions) == 0
                @constraint(model, x[v, s] == 0)
            end
        end
    end

    # every function should be processed (and only once)
   for s = 1:nS
        addition = 0
        for v = 1:nV
            addition += x[v, s]
        end
        @constraint(model, addition == 1)
    end

    # x[v, k, s] and x[v, k]
    for s = 1:nS
        for v = 1:nV
            @constraint(model, x[v, s] <= x0[v])
        end
    end

    # x[v,k] and y[v]
    for v = 1:nV
        @constraint(model, x0[v] <= y[v])
    end

    # node capacity constraint
    for v = 1:nV
            addition = 0
            for s = 1:nS
                addition += DataRate[kth] * x[v, s]
            end
            @constraint(model, addition == R0[v])
        @constraint(model, R0[v] <= Node_cap[v] * y[v])
    end

        # node capacity constraint
        # for v = 1:nV
        #     addition = 0
        #     for k = 1:nK
        #         for s = 1:nS
        #             addition += DataRate[k] * x[v, k, s]
        #         end
        #     end
        #     @constraint(model, addition <= Node_cap[v] * y[v])
        # end

    return model
end

function TrafficroutingMILPforsp(model, x, r, z, V, L, PairC, nV, nI, nL, nK, nS, nP, z0, kth)

    ## SFC constraints
        for s = 1:(nS+1)
            for i = 1:nI
                outgoinglinks = findall(L[:, 1] .== i)
                incominglinks = findall(L[:, 2] .== i)
                # addition2 = zeros(1, 1, 1, nP)
                # addition4 = zeros(1, 1, 1, nP)

                # for t = 1:length(outgoinglinks)
                #     addition2 .= addition2 .+ z[outgoinglinks[t], k, s, :]
                # end
                
                # for t = 1:length(incominglinks)
                #     addition4 .= addition4 .+ z[incominglinks[t], k, s, :]
                # end

                addition2 = sum(z[outgoinglinks, s, :], dims=1)
                addition4 = sum(z[incominglinks, s, :], dims=1)

                # println(size(addition2))
                
                v = sum((1:nV) .* (i .== V))
                
                # @constraint(model, addition2 .<= ones(1, 1, 1, nP))
                @constraint(model, addition2 .<= 1)
                
                if s == 1
                    if i == PairC[kth, 1]
                        # @constraint(model, addition2 .== ones(1, 1, 1, nP))
                        # @constraint(model, addition4 .== zeros(1, 1, 1, nP))
                        @constraint(model, addition2 .== 1)     # 情况1
                        @constraint(model, addition4 .== 0)
                    elseif v > 0
                        # @constraint(model, addition4 .- addition2 .== ones(1, 1, 1, nP) .* x[v, k, s])
                        @constraint(model, addition4 .- addition2 .== x[v, s])    # 情况2
                    else
                        @constraint(model, addition2 .== addition4)    # 情况6
                    end
                elseif s == (nS+1)
                    if i == PairC[kth, 2]
                        # @constraint(model, addition2 .== zeros(1, 1, 1, nP))
                        # @constraint(model, addition4 .== ones(1, 1, 1, nP))
                        @constraint(model, addition2 .== 0)                       # 情况5
                        @constraint(model, addition4 .== 1)
                    elseif v > 0
                        #@constraint(model, addition4 .- addition2 .== ones(1, 1, 1, nP) .* (-x[v, k, s-1]))
                        @constraint(model, addition4 .- addition2 .== (-x[v, s-1]))          # 情况4
                    else
                        @constraint(model, addition2 .== addition4)                # 情况6
                    end
                else
                    if v > 0
                        #@constraint(model, addition4 .- addition2 .== ones(1, 1, 1, nP) .* (x[v, k, s] - x[v, k, s-1]))
                        @constraint(model, addition4 .- addition2 .== (x[v, s] - x[v, s-1]))       # 情况3
                    else
                        @constraint(model, addition2 .== addition4)            # 情况6
                    end
                end
            end
        end

    # rz relation constraint
    
    @constraint(model, r[:, :, :] .<= z[:, :, :])
    

    
    for s = 1:(nS + 1)
        for i = 1:nI
            outgoinglinks = findall(L[:, 1] .== i)
            incominglinks = findall(L[:, 2] .== i)
            addition1 = 0
            addition2 = zeros(1, 1, nP)
            addition3 = 0
            addition4 = zeros(1, 1, nP)
            v = sum((1:nV) .* (i .== V))

            for t = 1:length(outgoinglinks)
                for p = 1:nP
                    addition1 += r[outgoinglinks[t], s, p]
                end
                addition2 = addition2 .+ r[outgoinglinks[t], s, :]
            end

            for t = 1:length(incominglinks)
                for p = 1:nP
                    addition3 += r[incominglinks[t], s, p]
                end
                addition4 = addition4 .+ r[incominglinks[t], s, :]

            end

            if s == 1
                if i == PairC[kth, 1]
                    @constraint(model, addition3 - addition1 == -1)      # 15 情况 1
                elseif v > 0
                    @constraint(model, addition3 - addition1 == x[v, s])   # 15 情况 2
                    @constraint(model, addition4 .- addition2 .<= x[v, s])   # 17
                else
                    @constraint(model, addition2 .== addition4)               # 16
                end
            elseif s == (nS + 1)
                if i == PairC[kth, 2]
                        @constraint(model, addition3 - addition1 == 1)    # 15 情况 5
                elseif v > 0
                    @constraint(model, addition3 - addition1 == -x[v, s - 1])     # 15 情况 4
                    @constraint(model, addition4 .- addition2 .>= -x[v, s - 1])    # 18   
                else
                    @constraint(model, addition2 .== addition4)         # 16
                end
            else
                if v > 0
                    @constraint(model, addition3 - addition1 == (x[v, s] - x[v, s - 1]))         # 15 情况 3   
                    @constraint(model, addition4 .- addition2 .<= x[v, s])       # 17
                    @constraint(model, addition4 .- addition2 .>= -x[v, s - 1])      # 18
                else
                    @constraint(model, addition2 .== addition4)         # 16
                end
            end
        end
    end
    

    # Valid inequalities
    for s = 1:(nS + 1)
        for l = 1:nL    
            addition = 0
            for p =1:nP
                addition = addition + r[l, s, p]
            end
            @constraint(model, addition <= z0[l])
        end
            # for p = 1:nP
            #     #@constraint(model, z[:, :, s, p] .<= z0[:, :])
            #     addition .= addition .+ r[:, :, s, p]
            # end
            #@constraint(model, addition .<= z0[:, :])
    end

    return model
end

function Linkcapforsp(model, r, R, DataRate, Link_cap, nL, nK, kth)
    # r R related constraint
    for l = 1:nL
        @constraint(model, R[l] == sum(sum(DataRate[kth] * r[l, :, :])))
    end

    # Link Capacity Constraint
    for l = 1:nL
        @constraint(model, R[l] <= Link_cap[l])
    end

    return model
end

function E2EreliabilityConsforsp(model, x, x0, r, z, z0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, nP, kth)
    # z z0 relation
    for s = 1:(nS + 1)
        #addition = zeros(nL, nK)
        for p = 1:nP
            @constraint(model, z[:, s, p] .<= z0[:])
           # addition .= addition .+ r[:, :, s, p]
        end
        # @constraint(model, addition .<= z0[:, :])
    end

    
        addition = 0
        for l = 1:nL
            addition = addition + log(Link_reliability[l]) * z0[l]
        end
        @constraint(model, addition .== rel_L)
        addition = 0
        for v = 1:nV
            addition = addition + log(Node_reliability[v]) * x0[v]
        end
        @constraint(model, addition .== rel_N)

        @constraint(model, rel_L + rel_N >= log(E2Ereliability[kth]))
    

    return model
end

function E2EdelayConsforsp(model, x, r, z, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, nP, kth)

        for s = 1:(nS + 1)
            addition2 = 0
            for p = 1:nP
                addition = 0
                for l = 1:nL
                    addition += Link_delay[l] * z[l, s, p]
                    addition2 += Link_delay[l] * r[l, s, p]
                end
                @constraint(model, addition <= theta[s])
            end
            @constraint(model, addition2 <= theta[s])
        end

        @constraint(model, sum(theta[:]) == theta_L)



        addition = 0
        for v = 1:nV
            addition += sum(x[v, :] .* Node_delay[v])
        end
        @constraint(model, addition == theta_N)
    

    
        @constraint(model, theta_L + theta_N <= E2EDelay[kth])
    

    return model
end

######################### 以下函数用于列生成算法子问题的Novel LP
function VNFplacementforspNovelLP(model, x, x0, y, V_processfunction, Profun, FunctionChainC, Node_cap, DataRate, nK, nS, nV, R0, kth)

    ## this function add constraints (1)-(4) in paper

    # no provide constraint (which doesn't appear in the paper) 
    for s = 1:nS
        for v = 2:nV      # here we always assume that node V[1] can process all service functions
            curprocessfunctions = V_processfunction[1:Profun, v]
            if sum(FunctionChainC[s, kth] .== curprocessfunctions) == 0
                @constraint(model, x[v, s] == 0)
            end
        end
    end

    # every function should be processed (and only once)
   for s = 1:nS
        addition = 0
        for v = 1:nV
            addition += x[v, s]
        end
        @constraint(model, addition == 1)
    end

    # x[v, k, s] and x[v, k]
    for s = 1:nS
        for v = 1:nV
            @constraint(model, x[v, s] <= x0[v])
        end
    end

    # x[v,k] and y[v]
    for v = 1:nV
        @constraint(model, x0[v] <= y[v])
    end

    # node capacity constraint
    for v = 1:nV
            addition = 0
            for s = 1:nS
                addition += DataRate[kth] * x[v, s]
            end
            @constraint(model, addition == R0[v])
        @constraint(model, R0[v] <= Node_cap[v] * y[v])
    end

        # node capacity constraint
        # for v = 1:nV
        #     addition = 0
        #     for k = 1:nK
        #         for s = 1:nS
        #             addition += DataRate[k] * x[v, k, s]
        #         end
        #     end
        #     @constraint(model, addition <= Node_cap[v] * y[v])
        # end

    return model
end

function TrafficroutingMILPforspNovelLP(model, x, r, V, L, PairC, nV, nI, nL, nK, nS, kth)

    ## SFC constraints
    
    for s = 1:(nS + 1)
        for i = 1:nI
            outgoinglinks = findall(L[:, 1] .== i)
            incominglinks = findall(L[:, 2] .== i)
            
            addition2 = 0
            
            addition4 = 0

            v = sum((1:nV) .* (i .== V))

            for t = 1:length(outgoinglinks)

                addition2 = addition2 + r[outgoinglinks[t], s]

            end

            for t = 1:length(incominglinks)

                addition4 = addition4 + r[incominglinks[t], s]

            end

            if s == 1
                if i == PairC[kth, 1]
                    @constraint(model, addition4 - addition2 == -1)      # 15 情况 1
                elseif v > 0
                    @constraint(model, addition4 - addition2 == x[v, s])   # 15 情况 2
                else
                    @constraint(model, addition2 == addition4)         # 15 情况 6
                end
            elseif s == (nS + 1)
                if i == PairC[kth, 2]
                        @constraint(model, addition4 - addition2 == 1)    # 15 情况 5
                elseif v > 0
                    @constraint(model, addition4 - addition2 == -x[v, s - 1])     # 15 情况 4
                else
                    @constraint(model, addition2 == addition4)         # 15 情况 6
                end
            else
                if v > 0
                    @constraint(model, addition4 - addition2 == (x[v, s] - x[v, s - 1]))         # 15 情况 3   
                else
                    @constraint(model, addition2 == addition4)         # 15 情况 6
                end
            end
        end
    end
    

    return model
end

function LinkcapforspNovelLP(model, r, R, DataRate, Link_cap, nL, nK, kth)
    # r R related constraint
    for l = 1:nL
        @constraint(model, R[l] == sum(DataRate[kth] * r[l, :]))
    end

    # Link Capacity Constraint
    for l = 1:nL
        @constraint(model, R[l] <= Link_cap[l])
    end

    return model
end

function E2EreliabilityConsforspNovelLP(model, x, x0, r, r0, rel_L, rel_N, Link_reliability, Node_reliability, E2Ereliability, nV, nL, nK, nS, kth)
    # z z0 relation
    for s = 1:(nS + 1)
        
        @constraint(model, r[:, s] .<= r0[:])

    end

    
        addition = 0
        for l = 1:nL
            addition = addition + log(Link_reliability[l]) * r0[l]
        end
        @constraint(model, addition .== rel_L)
        addition = 0
        for v = 1:nV
            addition = addition + log(Node_reliability[v]) * x0[v]
        end
        @constraint(model, addition .== rel_N)

        @constraint(model, rel_L + rel_N >= log(E2Ereliability[kth]))
    
    return model
end

function E2EdelayConsforspNovelLP(model, x, r, theta, theta_L, theta_N, Link_delay, Node_delay, E2EDelay, nV, nL, nK, nS, kth)

        for s = 1:(nS + 1)
            addition2 = 0

                for l = 1:nL
                    
                    addition2 += Link_delay[l] * r[l, s]
                end

            @constraint(model, addition2 == theta[s])
        end

        @constraint(model, sum(theta[:]) == theta_L)



        addition = 0
        for v = 1:nV
            addition += sum(x[v, :] .* Node_delay[v])
        end
        @constraint(model, addition == theta_N)
    

    
        @constraint(model, theta_L + theta_N <= E2EDelay[kth])
    

    return model
end

###############以下函数用于列生成算法主问题
function ConsforMP(model, t, y, c, dict_dict, Link_cap, Node_cap, nK, nL)
    
    consalpha = [ @constraint(model, sum(t[k, :]) == 1) for k in 1:nK ]
    
    Lhs1 = [zero(AffExpr) for _ in 1:nL]
    for l = 1:nL
        for k = 1: nK 
            for C = 1: c[k]
                #println(dict_dict[Symbol(k, C)][:R][l])
                #println("l = $l, k = $k, C = $C")
                # 输出字典的键和值
                #println("Keys of the dict dictionary: $(keys(dict_dict))")
                #println("????????????????????l$l,k$k,c$C")
                #println("Values of the solution dictionary: $(values(dict_dict))")
                Lhs1[l] = Lhs1[l] + dict_dict[(k, C)][:R][l] * t[k, C]
                #println(Lhs[l])
            end
        end 
    end

    consbeta = [ @constraint(model, Lhs1[l] <= Link_cap[l]) for l in 1:nL ]

    Lhs = [zero(AffExpr) for _ in 1:nV, _ in 1:nK]
    for v = 1:nV
        for k = 1:nK
            for C = 1:c[k]
                Lhs[v, k] += dict_dict[(k, C)][:x0][v] * t[k, C]
            end
        end
    end
    
    conspi = [ @constraint(model, Lhs[v, k] <= y[v]) for v in 1:nV, k in 1:nK]

    Lhs = [zero(AffExpr) for _ in 1:nV]
    for v = 1:nV
        for k = 1:nK
            for C = 1:c[k]
                Lhs[v] += dict_dict[(k, C)][:R0][v] * t[k, C]
                #print(Lhs[v])
            end
        end
    end 

    conseta = [ @constraint(model, Lhs[v] <= Node_cap[v] * y[v]) for v in 1:nV]

    return consalpha, consbeta, conspi, conseta, Lhs1, model
    
end

