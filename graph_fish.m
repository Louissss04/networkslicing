
clear
load('r_seed.mat');
rng(r_seed);

load data_case2.mat
% number of vertices 
nI=112;
% number of links
nL=440;
% length in the SFC
nS=4;
% number of cloud nodes
nV=6;
I = 1:nI;
% cloud nodes
V=[42, 51, 72, 77, 93, 96];
%links
L=zeros(nL,2);
%edge_cap
Link_cap=zeros(nL,1);
%edge_delay
Link_delay=zeros(nL,1);
capacityparam = 40;
% number of flows
nK=20;
%directory
Connectivity = 10;                                 
fish_network = ['fish' '-' num2str(nK) '-' num2str(capacityparam)]
%number of problems
nprob = 50;
%number of function types
nflow_types=4;
%FunctionChainC = zeros(nK,nS,nprob);
flow_combs=combntns([1:nflow_types],nS);
nflow_combs=combntns(nflow_types,nS);               % why
%PairC = zeros(nK,2,nprob);
V_processfunction=zeros(nflow_types,nV);

isLV = ones(nL, 1);   % 这里是为了把云节点v的相邻链路全激活 但逻辑错了
for l=1:nL
   for v=1:nV
      if ~(LinkInfo(l,1) == V(v) || LinkInfo(l,2) == V(v))
         isLV(l) = 0;
      end
   end
end


for t=1:nprob
   t
   DG=sparse(zeros(nI, nI));           % DG = delay of G
   RG=sparse(ones(nI, nI));            % RG = relieabilty of G  
    shortpathdistV=1e6;
    while (shortpathdistV == 1e6)
       G=sparse(zeros(nI, nI));
       isL = rand(nL,1) <= Connectivity/10;   % 激活百分之Connectivity的链路
       isL(nL-1,1) = 1;                       % 强制D(k)与唯一相邻的节点的双向链路激活
       isL(nL,1) = 1;
       isL = ((isLV + isL) > 0);
       tmpnL = sum(isL,1);                    % tmpnl计算了有多少链路被激活 本代码中应该就是440
       for l = 1:nL
          if isL(l) == 1                      % 如果链路被激活 G中存为1
              G(LinkInfo(l,1),LinkInfo(l,2)) = 1;
          end
       end
       for v=1:nV
         [tmpshortestpath,path,pred]=graphshortestpath(G,V(v),112);
         shortpathdistV = min(tmpshortestpath, shortpathdistV);
       end
    end
    tmpnL;
    originalgraphfile=['instances/fish/' fish_network '-' 'graph_info' '-' num2str(t)];
    fid=fopen(originalgraphfile,'w+');
    fprintf(fid,'%8d %8d %8d %8d\n',nI,tmpnL,nV,nflow_types);
    for l=1:nL
        if isL(l) == 1
            L(l,1) = LinkInfo(l,1);     % 这两行似乎没用
            L(l,2) = LinkInfo(l,2);
            Link_delay(l)=1+round(rand());
            Link_reliability(l)=0.994+0.001*round(5*rand()+0.5);
            Link_cap(l)= capacityparam*CapacityVector(l)/5;
        %%%%%%%%%%%%%%%%%%%%%%%%%%     need to look at this later 
            DG(LinkInfo(l,1),LinkInfo(l,2))=Link_delay(l);            % 两个权重图
            RG(LinkInfo(l,1),LinkInfo(l,2))=Link_reliability(l);
            fprintf(fid,'%8d %8d %8.4f %8.4f %8.4f\n',LinkInfo(l,1),LinkInfo(l,2),Link_cap(l),Link_delay(l),Link_reliability(l));
        end
    end
        
    Node_cap=200+ceil(400*rand(nV,1));
    Node_delay=3+round(3*rand(1,nV));
    Node_reliability=0.99+0.001*round(5*rand(nV,1)+0.5);
    for j=1:nV
       fprintf(fid,'%8d %8.2f %8.2f %8.2f \n',V(j),Node_cap(j),Node_delay(j),Node_reliability(j));
    end

    for j=1:nV
        tmp1=1:nflow_types;
        [tmp,index]=sort(rand(1,nflow_types));    % 生成一个随机数组并且排序 那么索引就是一个随机排列 用来模拟function chain的顺序
        V_processfunction(:,j) = tmp1(index);     
        V_nprocessfunction(j) = 2;
%        V_nprocessfunction(j) = 2 + round(rand());
    end
    for i=1:nflow_types
       for j=1:nV
          fprintf(fid,'%8d',V_processfunction(i,j));
       end
       fprintf(fid,'\n');
    end

    %save the graph data
    save(originalgraphfile,'nI','tmpnL','nV','V','V_processfunction','V_nprocessfunction',...
        'L','Link_delay','Link_cap',...
        'Link_reliability','Node_delay','Node_cap','Node_reliability');
    
    FunctionChainC = zeros(nK,nS);
    PairC = zeros(nK,2);

    
    % source destinction
    for j=1:nK
        tmp = round(0.5+(112-0.5)*rand());            % 随机生成一个备选起点
        shortpathdistV=1e6;
        while(sum(tmp==NoSource) > 0 || shortpathdistV == 1e6)
           tmp = round(0.5+(112-0.5)*rand());
           shortpathdistV=1e6;
           for v=1:nV
               [tmpshortestpath1,path,pred]=graphshortestpath(DG,tmp,V(v));     % why not G?
               [tmpshortestpath2,path,pred]=graphshortestpath(DG,V(v),112);
               shortpathdistV = min(max(tmpshortestpath1, tmpshortestpath2), shortpathdistV);  % 这一段保证了起点至少与一个云节点相连
           end
        end
        PairC(j,1)=tmp;            % PairC第一列存储了起点 第二列全部都是终点112
        PairC(j,2)=112;
    end
    
    % service function chain
    tmp=round(0.5+(nflow_combs-0.5)*rand(nK,1));   % 在这里tmp生成了1-nflowcombs等可能的随机数
    flowsfunctions=flow_combs(tmp',:);             % 先决定每个k选哪些function(无顺序)
    for j=1:nK                                     % 再随机打乱顺序存进FunctionChainC
        d=randperm(nS);
        tmp=flowsfunctions(j,:);
        FunctionChainC(j,:)=tmp(d);
    end
    
    
    DataRate=ceil(1+40*rand(nK,1));         % 对于同一个k, 所有的s速率都一样
    shortpathdist=zeros(nK,1);
    for i=1:nK
        [shortpathdist(i),path,pred]=graphshortestpath(DG,PairC(i,1),PairC(i,2));
        [Rshortpathdist(i),path,pred]=graphshortestpath(-log(RG),PairC(i,1),PairC(i,2));
    end
    E2Edelays = round(20+(3*shortpathdist)+5*rand());                    % 按照给定公式给出延迟和可靠性的阈值
    E2Ereliabilities=exp(-(4*Rshortpathdist-log(0.99)-log(0.99)));
    flowfile=['instances/fish/' fish_network '-' 'flow_info' '-' num2str(t)];
   
    fid=fopen(flowfile,'w+');
    fprintf(fid,'%8d %8d %8d \n',nK,nS,nflow_types);

    for i=1:nK
        fprintf(fid,'%12d %12d',PairC(i,1),PairC(i,2));
        fprintf(fid,'%12.4f',DataRate(i));
        fprintf(fid,'%12.4f',E2Edelays(i));
        fprintf(fid,'%12.4f',E2Ereliabilities(i));
        fprintf(fid,'\n');
    end

    for i=1:nS
        for j=1:nK
            fprintf(fid,'%4d',FunctionChainC(j,i));
        end
        fprintf(fid,'\n');
    end

    save(flowfile,'nK','nS','nflow_types','PairC','E2Edelays','E2Ereliabilities',...
        'DataRate','FunctionChainC');
    
  
    
end
















