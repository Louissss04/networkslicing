using FileIO
using Random


function readGraphFile(graphFile)

   # read the graph file
   io = readlines(graphFile)

   # read numbers: nI, nL, nV, nF
   str = split(io[1], (' ', '\t'), keepempty=false)
   nI = parse(Int, str[1])
   nL = parse(Int, str[2])
   nV = parse(Int, str[3])
   nF = parse(Int, str[4])
   #println(nI,"  ",nL,"  ",nV,"  ",nF)
   
   # read link information
   L = Array{Int64}(undef,nL,2);
   Link_cap = Array{Float64}(undef, nL)
   Link_delay = Array{Float64}(undef, nL)
   Link_reliability = Array{Float64}(undef, nL)
   l = 1;
   for line in io[2:(nL+1)]
      str = split(line, (' ', '\t'), keepempty=false)
      L[l,1] = parse(Int, str[1])
      L[l,2] = parse(Int, str[2])
      Link_cap[l] = parse(Float64, str[3])
      Link_delay[l] = parse(Float64, str[4])
      Link_reliability[l] = parse(Float64, str[5])
      l = l + 1
   end
   # println(L)
   # println(Link_cap) 
   # println(Link_delay)
   # println(Link_reliability)
 
   # read cloud node information
   V = Array{Int64}(undef, nV)
   Node_cap = Array{Float64}(undef, nV)
   Node_delay = Array{Float64}(undef, nV)
   Node_reliability = Array{Float64}(undef, nV)
   v = 1;
   for line in io[(nL+2):(nL+nV+1)]
      str = split(line, (' ', '\t'), keepempty=false)
      V[v] = parse(Int, str[1])
      Node_cap[v] = parse(Float64, str[2])
      Node_delay[v] = parse(Float64, str[3])
      Node_reliability[v] = parse(Float64, str[4])
      v = v + 1
   end
   # println(V)
   # println(Node_cap)
   # println(Node_delay)
   # println(Node_reliability)

   #read processing functions of cloud node
   V_processfunction = Array{Int64}(undef, nF, nV)
   f = 1;
   for line in io[(nL+nV+2):(nL+nV+nF+1)]
      str = split(line, (' ', '\t'), keepempty=false)
      for v in 1:nV
         V_processfunction[f, v] = parse(Int, str[v])
      end
      f = f + 1
   end
   # println(V_processfunction)

   outgoing_cap = zeros(nV)
   incoming_cap = zeros(nV)
   for v in 1:nV 
      outgoinglinks = findall(L[:,1].==V[v])
      incominglinks = findall(L[:,2].==V[v])
      outgoing_cap[v] = sum(Link_cap[outgoinglinks])
      incoming_cap[v] = sum(Link_cap[incominglinks])
   end

   return nI, nL, nV, nF, L, V, Link_cap, Link_reliability, Link_delay, Node_cap, Node_delay, Node_reliability, outgoing_cap, incoming_cap, V_processfunction
end


function readFlowFile(flowFile)
   #println("Flow file: ", flowFile)

   # read the graph file
   io = readlines(flowFile)

   # read numbers: nK, nS
   str = split(io[1], (' ', '\t'), keepempty=false)
   nK = parse(Int, str[1])
   nS = parse(Int, str[2])
   #println(nK,"  ",nS)

   # read source, destination, DataRates
   PairC = Array{Int64}(undef,nK,2);
   DataRate = Array{Float64}(undef,nK)
   E2EDelay = Array{Float64}(undef,nK)
   E2Ereliability = Array{Float64}(undef,nK)
   k = 1;
   for line in io[2:(nK+1)]
      str = split(line, (' ', '\t'), keepempty=false)
      PairC[k,1] = parse(Int, str[1])
      PairC[k,2] = parse(Int, str[2])
      DataRate[k] = parse(Float64, str[3])
      E2EDelay[k] = parse(Float64, str[4])
      E2Ereliability[k] = parse(Float64, str[5])
      k = k + 1
   end

    # read function chain information
   FunctionChainC = Array{Int64}(undef,nS,nK)
    s = 1
   for line in io[(nK+2):(nK+nS+1)]
   str = split(line, (' ', '\t'), keepempty=false)
      k = 1
      for k in 1:nK
         FunctionChainC[s, k] = parse(Int, str[k])
      end
      s = s + 1
   end
   return nK, nS, PairC, DataRate, FunctionChainC, E2EDelay, E2Ereliability
end


function readProb(graphFile, flowFile, ntests) 
 
   nI, nL, nV, nF, L, V, Link_cap, Link_reliability, Link_delay, Node_cap, Node_delay, Node_reliability, outgoing_cap, incoming_cap, V_processfunction = readGraphFile(graphFile)
   #println(nI, "\n", nL, "\n", nV, "\n", nF, "\n", L, "\n", V, "\n", 
   #        Link_cap, "\n", Node_cap, "\n", V_processfunction)

   nK, nS, PairC, DataRate, FunctionChainC, E2EDelay, E2Ereliability = readFlowFile(flowFile)
   #println(nK, "\n", nS, "\n", PairC, "\n", DataRate, "\n", FunctionChainC)
    
   # AssignCost = round.(20*rand(MersenneTwister(ntests),nV,nK,nS))
   # ActivationCost = round.(200*rand(MersenneTwister(ntests),nV))
   return nI, nL, nV, nF, L, V, Link_cap, Link_reliability, Link_delay, Node_cap, Node_delay, Node_reliability, outgoing_cap, incoming_cap, V_processfunction, nK, nS, PairC, DataRate, FunctionChainC, E2EDelay, E2Ereliability
end




