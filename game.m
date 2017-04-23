% Find a Minimum Vertex Cover through Game Theory
% This algorithm aims to verify an observation : 
% "When parameter R in the Snowdrift game satisfies R < 1/Kmax, all the 
% cooperators in a Nash equilibrium of a network game constitute a local
% minimum vertex cover."
% Chinese description of this observation :
% 当雪堆博弈的参数 R 满足 R < 1/Kmax 时,
% 网络博弈的纳什均衡中采用合作策略的节点构成极小节点覆盖.
% November 20, 2016, by HanzheTeng

clear variables
load graph.mat
n=length(graph);

% find the max degree of the input graph
% for symmetrical matrix, sum of rows equals to sum of columns
dmax=max(sum(graph)); 

% this observation is ture when r < 1/dmax in the snowdrift game
r=rand*(1/dmax);

% payoff matrix for the snowdrift game
payoff=[1,1-r;1+r,0];

%%%%% find a Nash equilibrium through evolutionary game %%%%%
% generate a stragety array randomly
% 1 = cooperator ; 2 = defector
stra=randi([1 2],1,n);
last_stra=zeros(1,n);
t=0;
while ~all(last_stra==stra)
    last_stra=stra;
    % compute the utility of each node
    utl=zeros(n,2); % (n,1)--utility when cooperate
                    % (n,2)--utility when defect
    for i=1:n
        for j=1:n
            if graph(i,j)==1
                utl(i,1)=utl(i,1)+payoff(1,stra(j)); % 1 = cooperate
                utl(i,2)=utl(i,2)+payoff(2,stra(j)); % 2 = defect
            end
        end
        % decide to cooperate or defect
        if utl(i,1) > utl(i,2) 
            stra(i)=1; 
        else
            stra(i)=2;
        end
    end
    t=t+1;
end
disp('1 = cooperator ; 2 = defector ');
disp(['strategy of each node : ' num2str(stra)]);
disp(['evolutionary times : ' num2str(t)]);

% observation: 
% When parameter R in the Snowdrift game satisfies R < 1/kmax, all the 
% cooperators in a Nash equilibrium of a network game constitute a local
% minimum vertex cover.

%%%%% verify this set is a local minimum vertex cover or not %%%%%
% cooperators constitute a vertex cover
index=find(stra==1); % indices of vertexes who are cooperators
set=graph;
set(index,:)=0;
set(:,index)=0;
if any(any(set))
    disp('Cooperators DON''T constitute a vertex cover.');
else
    disp('Cooperators constitute a vertex cover.');
end
% this vertex cover is local minimum
notcover=zeros(1,length(index));
m=1;
for i=index
    j=index(index~=i);
    set=graph;
    set(j,:)=0;
    set(:,j)=0;
    if any(any(set))
        notcover(m)=1;
    end
    m=m+1;
end
if all(notcover)
    disp('This vertex cover is local minimum.');
else
    disp('This vertex cover is NOT local minimum.');
end